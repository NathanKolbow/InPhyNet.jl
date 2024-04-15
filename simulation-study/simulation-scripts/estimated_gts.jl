# Full pipeline going from truenet > simulated gene trees > simulated sequences > inferred gene trees > inferred networks > netmerge

error("need to track the values we talk about in meeting w/ Kevin and Claudia and implement `save_estimated_gts_results`.")

if length(ARGS) != 6
    error("Usage: julia --project=X -tY estimated_gts.jl \"<true network abbreviation>\" <replicate number> <number of loci> <sequence length> <ils level (low/med/high)>")
end

###### Input parsing ######
netid = ARGS[1]
replicatenum = parse(Int64, ARGS[2])
ngt = parse(Int64, ARGS[3])
seq_len = parse(Int64, ARGS[4])
ils_level = ARGS[5]

(ngt == 100 || ngt == 1000) || error("Number of loci $ngt not allowed; must be 100 or 1,000.")
(seq_len == 100 || seq_len == 1000) || error("Sequence length $seq_len not allowed; must be 100 or 1,000.")
maxsubsetsize = 15
dmethod = "internode_count"
###########################
nhybrids = parse(Int64, split(netid, "r")[2])
data_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data"
###########################

include("helpers/helpers.jl")

# 1. gather ground truth network, constraint, distance matrix, and namelist
println("Loading data")
truenet = load_true_net_ils_adjusted(netid, replicatenum, ils_level)
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

# 8. InPhyNet inference
println("Merging networks")
mnet = nothing
try
    mnet = netnj(est_D, est_constraints, est_namelist)
catch e
end

# 9. Save results
# save_estimated_gts_results