include("helpers/pipelines.jl")

if length(ARGS) != 1
    error("Usage: julia perfect_subsets.jl \"<true network newick>\"")
end
truenewick = ARGS[1]

