function fixdir()
    currdir = split(pwd(), "/")
    while currdir[length(currdir)] != "network-merging"
        if length(currdir) <= 2
            cd(Base.source_dir())
        end
        cd("..")
        currdir = split(pwd(), "/")
    end
end
fixdir()

using Pkg; Pkg.activate(".")
cd("simulations")

using NetMerge, PhyloNetworks, StatsBase, Plots
include("plot-fxns.jl")

# DATA LOADING FUNCTIONS
function loadTrueData(netid::String, whichConstraints::Int64=1)
    if netid == "n40h4"
        truenet = readTopology("n40h4/n40h4.net")
        constraints = readMultiTopology("n40h4/true-constraints$(whichConstraints).net")
        return truenet, constraints
    else
        error("$(netid) not recognized")
    end
end


# PIPELINE FUNCTIONS
function runGroundTruthPipeline(truenet::HybridNetwork, constraints::Vector{HybridNetwork})
    # 1. Calculate distance matrix
    D, namelist = majorinternodedistance(truenet)

    # 2. Merge
    return netnj(D, constraints, namelist)
end


function runEstGtPipeline(estgts::Vector{HybridNetwork})
    # 1. MSCquartets
    error("MSCquartets not ported into pure julia yet")
    # hybsubsets, treesubset = ...

    # 2. Calculate distance matrix
    D, namelist = calculateAGID(estgts)

    # 3. SNaQ
    q, t = countquartetsintrees(estgts, showprogressbar=false)
    startingtree = PhyloNetworks.nj!(deepcopy(D), deepcopy(namelist))

    constraints = Array{HybridNetwork}(undef, length(hybsubsets))
    for (j, hybsub) in enumerate(hybusbsets)
        # Filter out all quartets that are not exclusively composed of taxa in `hybsub`
        temptaxonnumbers = [i for i in 1:length(t) if t[i] in hybsub]
        tempq = view(q, [i for i in 1:length(q) if all([number in temptaxonnumbers for number in q[i].taxonnumber])])
        tempdf = readTableCF(writeTableCF(tempq, t))

        # Estimate the network
        constraints[j] = snaq!(startingtree, tempdf, hmax=Int64(ceil(length(hybsub) / 3)), runs=10)
    end

    # 4. Merge
    return netnj(D, constraints, namelist=namelist)
end