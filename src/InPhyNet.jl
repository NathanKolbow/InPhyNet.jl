module InPhyNet

    include("imports.jl")
    include("constants.jl")
    include("Exceptions.jl")

    include("GraphHelpers.jl")

    include("inphynet/structs/SubNet.jl")
    include("inphynet/structs/ReticMap.jl")
    
    include("inphynet/compatibility.jl")
    include("inphynet/inphynet.jl")
    include("inphynet/internodedistance.jl")
    include("inphynet/mscquartetsinterface.jl")
    include("inphynet/subsetdecomp.jl")
    include("inphynet/subsetdecomp/njhierarchdecomp.jl")
    include("inphynet/subsetdecomp/satedecomp.jl")

    export inphynet,
        decomposeFromQuartets,
        majorinternodedistance, internodedistance, calculateAGID,
        majorinternodecount, internodecount, calculateAGIC,
        parsequartets, SubNet,
        findvalidpairs, findsiblingpairs, findoptQidx, ReticMap, updateconstraints!, Edge, mergeconstraintnodes!,
        njHierarchDecomp, pruneTruthFromDecomp,
        sateIdecomp, sateIIdecomp,
        SolutionDNEError, ConstraintError,
        are_compatible_heuristic, are_compatible_after_merge

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