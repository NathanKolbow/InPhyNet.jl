include("imports.jl")
include("constants.jl")
include("structs/InterimParentalTree.jl")
include("structs/Lineage.jl")
include("structs/LineageNode.jl")
include("algos/proba.jl")
include("algos/coalescent.jl")
include("algos/arbitrary-proba.jl")
include("algos/ParentalTrees.jl")
include("algos/ManipHelpers.jl")

export Lineage,
    LineageNode,
    InterimParentalTree, IPT,
    nlineages,
    lineages,
    getParentalTrees,
    ntaxa,
    coalesce