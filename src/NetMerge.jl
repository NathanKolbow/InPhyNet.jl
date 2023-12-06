module NetMerge

    include("imports.jl")
    include("constants.jl")

    include("GraphHelpers.jl")

    include("netmerge/structs/SubNet.jl")
    include("netmerge/structs/ReticMap.jl")
    
    include("netmerge/netnjmerge.jl")
    include("netmerge/internodedistance.jl")
    include("netmerge/mscquartetsinterface.jl")
    include("netmerge/subsetdecomp.jl")
    include("netmerge/subsetdecomp/nj-decomp.jl")

    export netnj, netnj!,
        decomposeFromQuartets,
        majorinternodedistance, internodedistance, calculateAGID,
        parsequartets, SubNet,
        findvalidpairs, findoptQidx, ReticMap, updateconstraints!, Edge, mergeconstraintnodes!

# include("ptrees/structs/InterimParentalTree.jl")
# include("ptrees/structs/Lineage.jl")
# include("ptrees/structs/LineageNode.jl")

# include("ptrees/proba.jl")
# include("ptrees/coalescent.jl")
# include("ptrees/arbitrary-proba.jl")
# include("ptrees/ParentalTrees.jl")
# include("ptrees/ManipHelpers.jl")

# export Lineage,
#     LineageNode,
#     InterimParentalTree, IPT,
#     nlineages,
#     lineages,
#     getParentalTrees,
#     ntaxa,
#     coalesce

end