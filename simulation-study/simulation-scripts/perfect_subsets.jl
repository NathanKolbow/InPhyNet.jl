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
if perfect_sims_already_performed(netid, replicatenum, maxsubsetsize)
    @info "Simulations already performed for $(netid)-$(replicatenum) w/ subset size $(maxsubsetsize); skipping."
    exit()
end

# 1. gather ground truth network, constraint, distance matrix, and namelist
truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, dmethod)

seed = parse(Int64, "$(truenet.numTaxa)42$(truenet.numHybrids)42$(replicatenum)")
Random.seed!(seed)

# 2. run robustness testing
println("- Running robustness testing for $(netid) ($(replicatenum)), max: $(maxsubsetsize)")
esterrors, majortreeRFs, gausserrors, constraintdiffs, nretics_est =
    monophyleticRobustness(truenet, constraints, D, namelist, nsim=nsim, displayprogress=true)
constraintdiffs = sum(constraintdiffs, dims=1)[1,:]

# 3. save results
savePerfectResults(
    truenet,
    constraints,
    esterrors,
    majortreeRFs,
    gausserrors,
    constraintdiffs,
    nretics_est,
    replicatenum,
    maxsubsetsize
)

#############
# Profiling #
#############

# include("helpers/helpers.jl")
# truenet, constraints, D, namelist = loadPerfectData("n100r5", 1, 15, "internode_count")
# _ = monophyleticRobustness(truenet, constraints, D, namelist, nsim=15)
# # Original:                 13.85   seconds
# # `findfirst` -> `nodemap`: 3.84    seconds

# using StatProfilerHTML
# @profilehtml esterrors, gausserrors, constraintdiffs, nretics_est = monophyleticRobustness(truenet, constraints, D, namelist, nsim=100)