using PhyloNetworks

# FUNCTIONS FOR LOADING THE DATA EXAMPLES USED IN THE DOCUMENTATION WALKTHROUGH
function load_inphynet_example_gts()::Vector{HybridNetwork}

    path = joinpath(dirname(@__FILE__), "..", "examples", "n50_example_gts.treefile")
    return readmultinewick(path)

end