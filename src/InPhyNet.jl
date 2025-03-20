
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
    include("inphynet/subsetdecomp/pruning.jl")
    include("inphynet/subsetdecomp/satedecomp.jl")

    include("inphynet/pipeline.jl")
    include("data_examples.jl")

    export inphynet,
        majorinternodedistance, internodedistance, calculateAGID,
        majorinternodecount, internodecount, calculateAGIC,
        updateconstraints!, Edge, mergeconstraintnodes!,
        prune_network, prune_networks,
        sateIdecomp, sateIIdecomp,
        centroid_edge_decomposition,
        SolutionDNEError, ConstraintError,
        are_compatible_heuristic, are_compatible_after_merge,
        # DOCS WALKTHROUGH DATA LOADING FXNS
        load_inphynet_example_gts

end