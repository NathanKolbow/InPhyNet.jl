include("../../helpers/helpers.jl")
include("try_all_roots_helpers.jl")

ENV["JULIA_DEBUG"] = Main

truenet, constraints, D, namelist = loadPerfectData("n50r2", 1, 15, "AGIC")
Random.seed!(42)
addnoise!(D, Normal(9., 9.))

netnj(D, constraints[1:2], namelist)


try_all_roots_netnj!(D, constraints[1:2], namelist)