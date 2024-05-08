# MAKE SURE TO RUN WITH `julia --project -tX ...`
# if running from a screen session: run from simulation-study/simulation-scripts/
# - example params: `julia --project=../.. -t4 ./perfect_subsets.jl n50r2 1 15 AGIC false false 100`

if length(ARGS) != 6 && length(ARGS) != 7
    error("Usage: julia perfect_subsets.jl <true network abbreviation> <replicate number> <maximum subset size> <distance method> <true/false (all have outgroup)> <true/false (outgroup removed after re-root)> [number of sims]")
end

# distance method options:
# - "internode_count"

###### Input parsing ######
netid = ARGS[1]
replicatenum = parse(Int64, ARGS[2])
maxsubsetsize = parse(Int64, ARGS[3])
dmethod = ARGS[4]
nsim = 1000
all_have_outgroup = parse(Bool, ARGS[5])
outgroup_removed_after_reroot = parse(Bool, ARGS[6])
if length(ARGS) == 7 nsim = parse(Int64, ARGS[7]) end
###########################
# netid, replicatenum, maxsubsetsize, dmethod, nsim = ("n200r10", 1, 25, "internode_count", 1000)

include("helpers/helpers.jl")
InPhyNet.TIEWARNING = true  # disables the warning message when there are ties

# 0. if the sims have already been performed, skip this
n_already_performed = n_perfect_sims_already_performed(netid, replicatenum, maxsubsetsize, all_have_outgroup, outgroup_removed_after_reroot)
nsim += 1 - n_already_performed     # +1 for the 0 error sims
if nsim <= 0
    @info "Simulations already performed for $(netid)-$(replicatenum) w/ subset size $(maxsubsetsize); skipping."
    exit()
end

# 1. gather ground truth network, constraint, distance matrix, and namelist
truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, dmethod, all_have_outgroup=all_have_outgroup, outgroup_removed_after_reroot=outgroup_removed_after_reroot)

seed = parse(Int64, "$(truenet.numTaxa)42$(truenet.numHybrids)42$(replicatenum)")
Random.seed!(seed)

# 2. run robustness testing
println("- Running robustness testing for $(netid) ($(replicatenum)), max: $(maxsubsetsize)")
esterrors, esterrors_without_missing_retics, majortreeRFs, gausserrors, constraintdiffs, nretics_est =
    monophyleticRobustness(truenet, constraints, D, namelist, nsim=nsim, displayprogress=false, do_no_noise_sim=(n_already_performed == 0))
constraintdiffs = sum(constraintdiffs, dims=1)[1,:]

# 3. filter out the results that had constraint errors
keep_idxs = esterrors .!= -2.

# 4. save results
savePerfectResults(
    truenet,
    constraints,
    esterrors[keep_idxs],
    esterrors_without_missing_retics[keep_idxs],
    majortreeRFs[keep_idxs],
    gausserrors[keep_idxs],
    constraintdiffs[keep_idxs],
    nretics_est[keep_idxs],
    replicatenum,
    maxsubsetsize,
    all_have_outgroup,
    outgroup_removed_after_reroot
)