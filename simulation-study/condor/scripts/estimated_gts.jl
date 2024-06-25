# Full pipeline going from true_net > simulated gene trees > simulated sequences > inferred gene trees > inferred networks > netmerge

###### Input parsing ######
netid = ARGS[1]
replicatenum = parse(Int64, ARGS[2])
ngt = parse(Int64, ARGS[3])
seq_len = parse(Int64, ARGS[4])
ils_level = ARGS[5]
maxsubsetsize = parse(Int64, ARGS[6])
nruns = 10

dmethod = "AGIC"
###########################
nhybrids = parse(Int64, split(netid, "r")[2])
data_dir = joinpath(pwd(), "data")
if !isdir(data_dir) mkdir(data_dir) end
if !isdir(joinpath(data_dir, "temp_data")) mkdir(joinpath(data_dir, "temp_data")) end
###########################

using Distributed
println("Loading helpers.jl")
# IMPORTANT: run `include` once on the main worker so that all the precompilation is done, then run `include`
#            with @everywhere so that each worker has access to necessary fxns
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/estimated_gts-helpers.jl")
@everywhere include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")
@everywhere include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/estimated_gts-helpers.jl")

# 0. check if we've already run these sims. if so, don't bother running again
if estimated_sims_already_performed(netid, replicatenum, ngt, seq_len, ils_level, maxsubsetsize)
    println("Simulations already performed for $(netid)-$(replicatenum) w/ ngt $(ngt), seq_len $(seq_len), ils_level $(ils_level); skipping.")
    exit()
else
    println("Parameter set: ($(netid), rep=$(replicatenum), ngt=$(ngt), len=$(seq_len), ils=$(ils_level), m=$(maxsubsetsize), d=$(dmethod))")
end

# 1/2. print some basic info
print_basic_info()

# 1. gather ground truth network, constraint, distance matrix, and namelist
log("DATA", "Loading ground truth data.")
true_net = load_true_net_ils_adjusted(netid, replicatenum, ils_level)
for e in true_net.edge
    if e.length == -1. e.length = 0.473 end
end

##############
# FILE PATHS #
##############
checkpoint_dir = joinpath(pwd(), "checkpoint_files")
if !isdir(checkpoint_dir) mkdir(checkpoint_dir) end

truegt_file = joinpath(checkpoint_dir, "truegt_$(netid)_$(replicatenum)_$(ngt)_$(ils_level).treefile")
seq_file_prefix = joinpath(checkpoint_dir, "seqfile_$(netid)_$(replicatenum)_$(ngt)_$(seq_len)_$(ils_level).phy")
estgt_file = joinpath(checkpoint_dir, "estgt_$(netid)_$(replicatenum)_$(ngt)_$(seq_len)_$(ils_level).treefile")
net_file = joinpath(checkpoint_dir, "estnets_$(netid)_$(replicatenum)_$(ngt)_$(seq_len)_$(ils_level)_$(maxsubsetsize)_$(dmethod).netfile")

if !isfile(estgt_file) touch(estgt_file) end
##############
seed = parse(Int64, "$(true_net.numTaxa)42$(true_net.numHybrids)42$(replicatenum)")



# 2. simulate gene trees (fast, do all at once)
true_gts::Vector{HybridNetwork} = simulate_gene_trees(truegt_file, ngt, seed)

# 3. true gt --> seq --> est gt (slow, do one at a time and use checkpointing)
pmap(
    (i, gt) -> est_gt_from_true_gt(gt, "$(seq_file_prefix)_$(i)", "$(estgt_file)_$(i)", data_dir, i),
    1:length(true_gts), true_gts
)
flush(stdout)
est_gts = [readTopology("$(estgt_file)_$(i)") for i=1:length(true_gts)]

# 4. estimate an NJ tree
est_D, est_namelist, nj_tre = estimate_nj_tree(est_gts)

# 5. decompose into disjoint subsets
subsets = subset_decomp(nj_tre, maxsubsetsize)

# 6. infer constraints w/ SNaQ
est_constraints, est_constraint_runtimes = snaq_constraints(est_gts, net_file, subsets, true_net, nj_tre, nruns, seed)

# 7. InPhyNet inference
log("InPhyNet", "Merging networks.")
try
    inphynet_time = @elapsed mnet = netnj(est_D, est_constraints, est_namelist)
    log("InPhyNet", "Merge successful.", :green)

    log("InPhyNet", "Saving results.")
    save_estimated_gts_results(netid, true_net, replicatenum, ngt,
        ils_level, maxsubsetsize, dmethod, seq_len, mnet, est_constraints,
        est_gts, est_constraint_runtimes, inphynet_time)
catch e
    log("InPhyNet", "Failed ($(typeof(e))).", :red)
end

# 8. print all the necessary info into the output file in case it doesn't get saved properly
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

