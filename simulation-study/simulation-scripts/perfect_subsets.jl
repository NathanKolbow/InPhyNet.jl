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
###########################

include("helpers/helpers.jl")

# 1. gather ground truth network, constraint, distance matrix, and namelist
truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, dmethod)

# 2. run robustness testing
esterrors, gausserrors, constraintdiffs, nretics_est = monophyleticRobustness(truenet, constraints, D, namelist, nsim=nsim)
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