# MAKE SURE TO RUN WITH `julia --project -tX ...`
# if running from a screen session: run from simulation-study/simulation-scripts/
# - `julia --project=../.. -tX ./perfect_subsets.jl ...`

#error("Double check in Slack how we're choosing the Gaussian standard error and refactor `monophyleticRobustness` to use this schema.")
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

@everywhere netid = $netid
@everywhere replicatenum = $replicatenum
@everywhere maxsubsetsize = $maxsubsetsize
@everywhere dmethod = $dmethod
@everywhere nsim = $nsim
###########################


# 1. gather ground truth network, constraint, distance matrix, and namelist
using Distributed
@everywhere begin
    include("helpers/helpers.jl")
    InPhyNet.TIEWARNING = true

    truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, dmethod)
end

# 2. run robustness testing
println("- Running robustness testing for $(netid) ($(replicatenum)), max: $(maxsubsetsize)")
esterrors, gausserrors, constraintdiffs, nretics_est =
    monophyleticRobustnessDistributed(truenet, constraints, D, namelist, nsim=nsim)
constraintdiffs = sum(constraintdiffs, dims=1)[1,:]

# 3. save results
savePerfectResults(
    truenet,
    constraints,
    esterrors,
    gausserrors,
    constraintdiffs,
    nretics_est,
    replicatenum
)


# Let's try using multiple CPUs...
using Distributed


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