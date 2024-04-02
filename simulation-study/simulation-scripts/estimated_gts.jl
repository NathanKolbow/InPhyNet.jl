# Full pipeline going from truenet > simulated gene trees > simulated sequences > inferred gene trees > inferred networks > netmerge

if length(ARGS) != 5
    error("Usage: julia estimated_gts.jl \"<true network abbreviation>\" <replicate number> <maximum subset size> \"<distance method>\", <number of gene trees>")
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

# 1.1 set seed
seed = parse(Int64, "$(truenet.numTaxa)42$(truenet.numHybrids)42$(replicatenum)")

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
seq_files = simulate_sequence_data(gts, truegt_file, seq_file)

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

# 7. infer constraints w/ SNaQ
println("Inferring networks")
net_file = joinpath(data_dir, "estnets_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).netfile")
Random.seed!(seed)
infer_constraints(estgt_file, net_file, subsets, nhybrids)

# 8. InPhyNet inference
println("Merging networks")
est_constraints = readMultiTopology(net_file)

mnet = netnj(est_D, est_constraints, est_namelist)