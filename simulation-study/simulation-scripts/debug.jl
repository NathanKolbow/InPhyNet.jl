# TEST RUNS
using InPhyNet, PhyloNetworks
include("helpers/helpers.jl")
include("helpers/simple-plot-fxns.jl")

# All done w/ replicate 1, maxsubsetsize 15, dmethod internode_count
# Format: netid (proportion drops from 1, proportion hits 0):
# - n50r2       (0.48, 0.92) - prop-1: 0.59
# - n50r5       (0.52, 1.07) - prop-1: 0.49
# - n100r5      (0.49, 0.91) - prop-1: 0.58
# - n100r10     (0.02, 0.78) - prop-1: 0.89
#   - repnum 4  (0.40, 0.84) - prop-1: 0.63
#   - repnum 6  (0.30, 0.71) - prop-1: 0.58
netid = "n100r10"
replicatenum = 1
maxsubsetsize = 15
dmethod = "internode_count"
nsim = 500

truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, dmethod)
esterrors, gausserrors, constraintdiffs, nretics_est =
    monophyleticRobustness(truenet, constraints, D, namelist, nsim=nsim, displayprogress=true)
constraintdiffs = sum(constraintdiffs, dims=1)[1,:]

plot_negative_one_prop(netid, upperTriangStd(D))

# DEBUG
using Revise, InPhyNet, PhyloNetworks
include("helpers/helpers.jl")

