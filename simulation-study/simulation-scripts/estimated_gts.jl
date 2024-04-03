# Full pipeline going from truenet > simulated gene trees > simulated sequences > inferred gene trees > inferred networks > netmerge

if length(ARGS) != 5
    error("Usage: julia --project=X -tY estimated_gts.jl \"<true network abbreviation>\" <replicate number> <maximum subset size> \"<distance method>\", <number of gene trees>")
end

###### Input parsing ######
netid = ARGS[1]
replicatenum = parse(Int64, ARGS[2])
maxsubsetsize = parse(Int64, ARGS[3])
dmethod = ARGS[4]
ngt = parse(Int64, ARGS[5])
###########################
nhybrids = parse(Int64, split(netid, "r")[2])
data_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data"
###########################

include("helpers/helpers.jl")

# 1. gather ground truth network, constraint, distance matrix, and namelist
println("Loading data")
truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, dmethod)
seed = parse(Int64, "$(truenet.numTaxa)42$(truenet.numHybrids)42$(replicatenum)")
for e in truenet.edge
    if e.length == -1. e.length = 0.473 end
end

# 2. simulate gene trees
println("Simulating gene trees")
using PhyloCoalSimulations
truegt_file = joinpath(data_dir, "truegt_$(netid)_$(replicatenum)_$(ngt).treefile")
gts = nothing

if isfile(truegt_file)
    gts = readMultiTopology(truegt_file)
else
    Random.seed!(seed)
    gts = simulatecoalescent(truenet, ngt, 1)
    writeMultiTopology(gts, truegt_file)
end

# 3. simulate sequences w/ seq-gen
println("Simulating sequences")
seq_file_prefix = joinpath(data_dir, "seqfile_$(netid)_$(replicatenum)_$(ngt).phy")
Random.seed!(seed)
seq_files = simulate_sequence_data(gts, truegt_file, seq_file_prefix)

# 4. infer gene trees w/ iqtree
println("Inferring gene trees")
estgt_file = joinpath(data_dir, "estgt_$(netid)_$(replicatenum)_$(ngt).treefile")
Random.seed!(seed)
estimate_gene_trees(seq_files, estgt_file)

# 5. estimate an NJ tree for subset decomposition
println("Estimating NJ tree")
est_gts = readMultiTopology(estgt_file)
est_D, est_namelist = calculateAGIC(est_gts)
nj_df = DataFrame(est_D, est_namelist)
nj_tre = nj(nj_df)

# 6. decompose into disjoint subsets
println("Decomposing taxa into disjoint subsets")
subsets = sateIdecomp(nj_tre, maxsubsetsize)
if minimum([length(s) for s in subsets]) < 4
    error("Smallest subset must have at least 4 taxa for SNaQ.")
end

# 7. infer constraints w/ SNaQ
println("Inferring networks")
net_file = joinpath(data_dir, "estnets_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).netfile")
Random.seed!(seed)
est_constraints = infer_constraints(estgt_file, net_file, subsets, nhybrids)

try
    # 8. InPhyNet inference
    println("Merging networks")

    mnet = netnj(est_D, est_constraints, est_namelist)

    # 9. Save results
    #### ----> WE NEED TO SAVE THE RUNTIME FOR EACH INDIVIDUAL SNAQ RUN SO WE CAN COMPUTE T_serial AND T_parallel
    ####       i.e. one column in the CSV will have entries looking like "[10, 20, 13, 5]"
    writeTopology(mnet, "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/data/output/estimated_gt_output/finalnet_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).netfile")
catch e
    open("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/data/output/estimated_gt_output/finalnet_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).netfile", "w+") do f
        write(f, "FAILED")
    end
end

writeMultiTopology(est_constraints, "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/data/output/estimated_gt_output/estconstraints_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).netfile")