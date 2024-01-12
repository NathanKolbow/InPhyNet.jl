# MAKE SURE TO RUN WITH `julia --project -tX ...`
# TODO: switch up how network newick is passed. passing as a string will probably get too long and cause problems eventuallyget
if length(ARGS) != 3
    error("Usage: julia perfect_subsets.jl \"<true network newick>\" <maximum subset size> \"<distance method>\"")
end
truenewick = ARGS[1]
maxsubsetsize = parse(Int64, ARGS[2])
dmethod = ARGS[3]

include("helpers/robustness-fxns.jl")
include("helpers/save-results.jl")

# 0. gather ground truth distance matrix & namelist
truenet = readTopology(truenewick)

D, namelist = (nothing, nothing)
if dmethod == "internode_count"
    D, namelist = majorinternodedistance(truenet)
else
    error("Unrecognized distance method specified.")
end

# 1. subset decomp w/ SATe-I decomp method (SATe-II doesn't seem to work very well for our purposes)
taxa_subsets = sateIdecomp(majorTree(truenet), maxsubsetsize)
constraints = pruneTruthFromDecomp(truenet, taxa_subsets)

# 2. run robustness testing
esterrors, gausserrors, constraintdiffs = monophyleticRobustness(truenet, constraints, D, namelist)
constraintdiffs = sum(constraintdiffs, dims=1)[1,:]

# 3. save results
