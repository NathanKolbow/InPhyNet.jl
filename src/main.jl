include("imports.jl")
include("constants.jl")
include("structs/Lineage.jl")
include("structs/LineageNode.jl")
include("algos/ParentalTrees.jl")
include("algos/ManipHelpers.jl")

export Lineage,
    LineageNode,
    nlineages,
    lineages,
    getParentalTrees,
    ntaxa,
    coalesce