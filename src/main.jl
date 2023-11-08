# module netmerging

include("imports.jl")
include("constants.jl")
include("structs/InterimParentalTree.jl")
include("structs/Lineage.jl")
include("structs/LineageNode.jl")
include("algos/ptrees/proba.jl")
include("algos/ptrees/coalescent.jl")
include("algos/ptrees/arbitrary-proba.jl")
include("algos/ptrees/ParentalTrees.jl")
include("algos/ptrees/ManipHelpers.jl")

export Lineage,
    LineageNode,
    InterimParentalTree, IPT,
    nlineages,
    lineages,
    getParentalTrees,
    ntaxa,
    coalesce

# end module