# MAKE SURE TO RUN WITH `julia --project -tX ...`
# if running from a screen session: run from simulation-study/simulation-scripts/
# - `julia --project=../.. -tX ./perfect_subsets.jl ...`
# - example params: `julia ... ./perfect_subsets.jl n50r2 1 15 internode_count 100`

if length(ARGS) != 4 && length(ARGS) != 5
    error("Usage: julia perfect_subsets.jl \"<true network abbreviation>\" <replicate number> <maximum subset size> \"<distance method>\" [number of sims]")
end

# distance method options:
# - "internode_count"

###### Input parsing ######
netid = ARGS[1]
replicatenum = parse(Int64, ARGS[2])
maxsubsetsize = parse(Int64, ARGS[3])
dmethod = ARGS[4]
nsim = 1000
if length(ARGS) == 5 nsim = parse(Int64, ARGS[5]) end
###########################
# netid, replicatenum, maxsubsetsize, dmethod, nsim = ("n200r10", 1, 25, "internode_count", 1000)

include("helpers/helpers.jl")
InPhyNet.TIEWARNING = true  # disables the warning message when there are ties

# 0. if the sims have already been performed, skip this
nsim += 1 - n_perfect_sims_already_performed(netid, replicatenum, maxsubsetsize)    # +1 for the 0 error sims
if nsim <= 0
    @info "Simulations already performed for $(netid)-$(replicatenum) w/ subset size $(maxsubsetsize); skipping."
    exit()
end

# 1. gather ground truth network, constraint, distance matrix, and namelist
truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, dmethod)

seed = parse(Int64, "$(truenet.numTaxa)42$(truenet.numHybrids)42$(replicatenum)")
Random.seed!(seed)

# 2. run robustness testing
println("- Running robustness testing for $(netid) ($(replicatenum)), max: $(maxsubsetsize)")
esterrors, esterrors_without_missing_retics, majortreeRFs, gausserrors, constraintdiffs, nretics_est =
    monophyleticRobustness(truenet, constraints, D, namelist, nsim=nsim, displayprogress=true)
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
    maxsubsetsize
)