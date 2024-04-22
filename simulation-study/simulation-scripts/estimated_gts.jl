# Full pipeline going from true_net > simulated gene trees > simulated sequences > inferred gene trees > inferred networks > netmerge


if length(ARGS) == 0
    # Inputs for testing
    @warn "No input arguments provided, running script w/ test parameters (netid = n200r10, replicatenum = 2, ngt = seq_len = 1000, ils_level = med)"
    push!(ARGS, "n200r10")
    push!(ARGS, "2")
    push!(ARGS, "1000")
    push!(ARGS, "1000")
    push!(ARGS, "med")
elseif length(ARGS) != 5
    error("Usage: julia --project=X -tY estimated_gts.jl \"<true network abbreviation>\" <replicate number> <number of loci> <sequence length> <ils level (low/med/high)>")
end

###### Input parsing ######
netid = ARGS[1]
replicatenum = parse(Int64, ARGS[2])
ngt = parse(Int64, ARGS[3])
seq_len = parse(Int64, ARGS[4])
ils_level = ARGS[5]

(ngt == 100 || ngt == 1000 || ngt == 5000) || error("Number of loci $ngt not allowed; must be 100 or 1,000.")
(seq_len == 500 || seq_len == 1000) || error("Sequence length $seq_len not allowed; must be 100 or 1,000.")
maxsubsetsize = 15
dmethod = "AGIC"
###########################
nhybrids = parse(Int64, split(netid, "r")[2])
data_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data"
###########################

include("helpers/helpers.jl")

# 0. check if we've already run these sims. if so, don't bother running again
if estimated_sims_already_performed(netid, replicatenum, ngt, seq_len, ils_level)
    @info "Simulations already performed for $(netid)-$(replicatenum) w/ ngt $(ngt), seq_len $(seq_len), ils_level $(ils_level); skipping."
    exit()
end

# 1. gather ground truth network, constraint, distance matrix, and namelist
@info "Loading data"
true_net = load_true_net_ils_adjusted(netid, replicatenum, ils_level)
seed = parse(Int64, "$(true_net.numTaxa)42$(true_net.numHybrids)42$(replicatenum)")
for e in true_net.edge
    if e.length == -1. e.length = 0.473 end
end

# 2. simulate gene trees
@info "Simulating gene trees"
using PhyloCoalSimulations
truegt_file = joinpath(data_dir, "truegt_$(netid)_$(replicatenum)_$(ngt).treefile")
gts = nothing

if isfile(truegt_file)
    gts = readMultiTopology(truegt_file)
else
    Random.seed!(seed)
    gts = simulatecoalescent(true_net, ngt, 1)
    writeMultiTopology(gts, truegt_file)
end

# 3. simulate sequences w/ seq-gen
@info "Simulating sequences"
seq_file_prefix = joinpath(data_dir, "seqfile_$(netid)_$(replicatenum)_$(ngt).phy")
Random.seed!(seed)
seq_files = simulate_sequence_data(gts, truegt_file, seq_file_prefix)

# 4. infer gene trees w/ iqtree
@info "Inferring gene trees"
estgt_file = joinpath(data_dir, "estgt_$(netid)_$(replicatenum)_$(ngt).treefile")
Random.seed!(seed)
estimate_gene_trees(seq_files, estgt_file)

# 5. estimate an NJ tree for subset decomposition
@info "Estimating NJ tree"
est_gts = readMultiTopology(estgt_file)
est_D, est_namelist = calculateAGIC(est_gts)
nj_df = DataFrame(est_D, est_namelist)
nj_tre = nj(nj_df)

# 6. decompose into disjoint subsets
@info "Decomposing taxa into disjoint subsets"
subsets = sateIdecomp(nj_tre, maxsubsetsize)
if minimum([length(s) for s in subsets]) < 4
    error("Smallest subset must have at least 4 taxa for SNaQ.")
end

# 7. infer constraints w/ SNaQ
@info "Inferring networks"
net_file = joinpath(data_dir, "estnets_$(netid)_$(replicatenum)_$(maxsubsetsize)_$(dmethod)_$(ngt).netfile")
Random.seed!(seed)
est_constraints, est_constraint_runtimes = infer_constraints(estgt_file, net_file, subsets, true_net)

# 8. InPhyNet inference
@info "Merging networks"
mnet = nothing
inphynet_time = @elapsed try
    global mnet
    mnet = netnj(est_D, est_constraints, est_namelist)
    @info "Network merge step succeeded"
catch e
    @info "Network merge step failed: $(typeof(e))"
end

# 9. Save results
@info "Saving results"
save_estimated_gts_results(netid, true_net, replicatenum, ngt,
    ils_level, maxsubsetsize, dmethod, seq_len,
    mnet, est_constraints, est_gts, est_constraint_runtimes, inphynet_time)