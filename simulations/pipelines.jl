function fixdir()
    currdir = splitdir(pwd())
    while currdir[length(currdir)] != "network-merging"
        if length(currdir) <= 2
            cd(Base.source_dir())
        end
        cd("..")
        currdir = splitdir(pwd())
    end
end
fixdir()

using Pkg; Pkg.activate(".")
cd("simulations")

using NetMerge, PhyloNetworks, StatsBase, Plots, StatsPlots, DataFrames, CSV
include("robustness-fxns.jl")
include("plot-fxns.jl")

# DATA LOADING FUNCTIONS
function loadTrueData(netid::String, whichConstraints::Int64=1)
    if netid == "n40h4"
        truenet = readTopology("n40h4/n40h4.net")
        constraints = readMultiTopology("n40h4/true-constraints$(whichConstraints).net")
        return truenet, constraints
    elseif netid == "n80h8"
        truenet = readTopology("n80h8/n80h8.net")
        constraints = readMultiTopology("n80h8/true-constraints$(whichConstraints).net")
        return truenet, constraints
    else
        error("$(netid) not recognized")
    end
end


# PIPELINE FUNCTIONS
function runGroundTruthRobustnessPipeline(truenet::HybridNetwork, constraints::Vector{HybridNetwork}=Vector{HybridNetwork}([]); nsim::Int64=100)
    if length(constraints) == 0
        decomp = njHierarchDecomp(majorTree(truenet), 16)
        constraints = pruneTruthFromDecomp(truenet, decomp)
    end

    # 1. Ground truth estimation
    D, namelist = majorinternodedistance(truenet)
    mnet = netnj(D, constraints, namelist)

    # 2. Distance matrix robustness tests
    robD1 = robustGauss(truenet, constraints, μ=1., σ=1., nsim=nsim)
    robD2 = robustGauss(truenet, constraints, μ=2., σ=2., nsim=nsim)
    robD3 = robustGauss(truenet, constraints, μ=4., σ=4., nsim=nsim)

    # 3. NNI robustness tests
    a0, b0, c0 = robustNNI(truenet, constraints, repeat([0.5], length(constraints)), nsim=2*nsim)
    a1, b1, c1 = robustNNI(truenet, constraints, repeat([1], length(constraints)), nsim=nsim)
    a2, b2, c2 = robustNNI(truenet, constraints, repeat([2], length(constraints)), nsim=nsim)
    a3, b3, c3 = robustNNI(truenet, constraints, repeat([3], length(constraints)), nsim=nsim)
    b0 = sum(b0, dims=1)[1,:]
    b1 = sum(b1, dims=1)[1,:]
    b2 = sum(b2, dims=1)[1,:]
    b3 = sum(b3, dims=1)[1,:]
    c0 = sum(c0, dims=1)[1,:]
    c1 = sum(c1, dims=1)[1,:]
    c2 = sum(c2, dims=1)[1,:]
    c3 = sum(c3, dims=1)[1,:]

    r1 = getNetDistances(truenet, mnet)
    r2 = DataFrame(
        gaussMean=[repeat([1.], length(robD1)); repeat([2.], length(robD2)); repeat([4.], length(robD3))],
        gaussStd=[repeat([1.], length(robD1)); repeat([2.], length(robD2)); repeat([4.], length(robD3))],
        dists=vcat(robD1, robD2, robD3)
    )
    r3 = DataFrame(
        nniMoves=[repeat([0.5], length(a0)); repeat([1], length(a1)); repeat([2], length(a2)); repeat([3], length(a3))],
        edgeheights=vcat(c0, c1, c2, c3),
        constraintdists=vcat(b0, b1, b2, b3),
        estdists=vcat(a0, a1, a2, a3)
    )

    return (
        r1,
        r2,
        r3
    )
end


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