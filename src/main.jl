include("structs/Lineage.jl")
include("structs/LineageNode.jl")
include("algos/ParentalTrees.jl")

export Lineage,
    LineageNode,
    nlineages,
    lineages,
    getParentalTrees,
    ntaxa,
    coalesce