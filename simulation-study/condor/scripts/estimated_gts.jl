# Full pipeline going from true_net > simulated gene trees > simulated sequences > inferred gene trees > inferred networks > netmerge

###### Input parsing ######
netid = ARGS[1]
replicatenum = parse(Int64, ARGS[2])
ngt = parse(Int64, ARGS[3])
seq_len = parse(Int64, ARGS[4])
ils_level = ARGS[5]
maxsubsetsize = parse(Int64, ARGS[6])

dmethod = "AGIC"
###########################
nhybrids = parse(Int64, split(netid, "r")[2])
data_dir = joinpath(pwd(), "data")
mkdir(data_dir)
mkdir(joinpath(data_dir, "temp_data"))
###########################

println("Loading helpers.jl")
# IMPORTANT: run `include` once on the main worker so that all the precompilation is done, then run `include`
#            with @everywhere so that each worker has access to necessary fxns
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")
@everywhere include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")

# 0. check if we've already run these sims. if so, don't bother running again
if estimated_sims_already_performed(netid, replicatenum, ngt, seq_len, ils_level, maxsubsetsize)
    println("Simulations already performed for $(netid)-$(replicatenum) w/ ngt $(ngt), seq_len $(seq_len), ils_level $(ils_level); skipping.")
    exit()
else
    println("Parameter set: ($(netid), rep=$(replicatenum), ngt=$(ngt), len=$(seq_len), ils=$(ils_level), m=$(maxsubsetsize), d=$(dmethod))")
end

# 1. gather ground truth network, constraint, distance matrix, and namelist
println("Loading data")
true_net = load_true_net_ils_adjusted(netid, replicatenum, ils_level)
seed = parse(Int64, "$(true_net.numTaxa)42$(true_net.numHybrids)42$(replicatenum)")
for e in true_net.edge
    if e.length == -1. e.length = 0.473 end
end

# 2. simulate gene trees
println("Simulating gene trees")
using PhyloCoalSimulations
truegt_file = joinpath(data_dir, "truegt_$(netid)_$(replicatenum)_$(ngt)_$(ils_level).treefile")
gts = nothing

if isfile(truegt_file)
    gts = readMultiTopology(truegt_file)
else
    Random.seed!(seed)
    gts = simulatecoalescent(true_net, ngt, 1)
    writeMultiTopology(gts, truegt_file)
end

# 3. simulate sequences w/ seq-gen
println("Simulating sequences")
seq_file_prefix = joinpath(data_dir, "seqfile_$(netid)_$(replicatenum)_$(ngt)_$(seq_len)_$(ils_level).phy")
estgt_file = joinpath(data_dir, "estgt_$(netid)_$(replicatenum)_$(ngt)_$(seq_len)_$(ils_level).treefile")
Random.seed!(seed)
seq_files = simulate_sequence_data(Vector{HybridNetwork}(gts), truegt_file, seq_file_prefix, estgt_file, data_dir)

# 4. infer gene trees w/ iqtree
println("Inferring gene trees")
Random.seed!(seed)
estimate_gene_trees(seq_files, estgt_file)

# 4.1 delete the simulated sequence data, that way we can save some space...
for seq_file in seq_files
    if isfile(seq_file)
        rm(seq_file)
    end
end

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
net_file = joinpath(data_dir, "estnets_$(netid)_$(replicatenum)_$(ngt)_$(seq_len)_$(ils_level)_$(maxsubsetsize)_$(dmethod).netfile")
Random.seed!(seed)
est_constraints, est_constraint_runtimes = infer_constraints(estgt_file, net_file, subsets, true_net)

# 8. InPhyNet inference
println("Merging networks")
mnet = nothing
try
    global mnet
    inphynet_time = @elapsed mnet = netnj(est_D, est_constraints, est_namelist)
    println("Network merge step succeeded")

    # 9. Save results
    println("Saving results")
    save_estimated_gts_results(netid, true_net, replicatenum, ngt,
        ils_level, maxsubsetsize, dmethod, seq_len,
        mnet, est_constraints, est_gts, est_constraint_runtimes, inphynet_time)
catch e
    println("Network merge step failed: $(typeof(e))")
end

# 10. print all the necessary info into the output file in case it doesn't get saved properly
println(netid)
println(writeTopology(true_net))
println(replicatenum)
println(ngt)
println(ils_level)
println(maxsubsetsize)
println(dmethod)
println(seq_len)
println(length(est_constraints))
for c in est_constraints
    println(writeTopology(c))
end
println(length(est_gts))
for t in est_gts
    println(writeTopology(t))
end
for r in est_constraint_runtimes
    println(r)
end