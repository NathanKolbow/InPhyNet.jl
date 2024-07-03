using InPhyNet, PhyloNetworks, StatsBase, DataFrames, CSV, Random
include("robustness-fxns.jl")

# FILE PATH HELPERS
function getBaseDir()
    basedir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/"
    if !isdir(basedir)
        basedir = "C:\\Users\\Nathan\\repos\\network-merging\\simulation-study\\"
    end
    if !isdir(basedir)
        basedir = "/Users/nkolbow/repos/network-merging/simulation-study/"
    end
    return basedir
end

getDataDir() = joinpath(getBaseDir(), "data")
getNetworkFilepath(netid::String) = joinpath(getDataDir(), "networks", "$(netid).netfile")
getNetworkFilepathLevel1(ntaxa::Int64) = joinpath(getDataDir(), "networks-level1", "n$(ntaxa).netfile")
getOutputFilepath(truenet::HybridNetwork) = joinpath(getDataDir(), "output", "n$((truenet.numTaxa)-1)r$(truenet.numHybrids).csv")
getOutputFilepathLevel1(truenet::HybridNetwork) = joinpath(getDataDir(), "networks-level1", "n$(truenet.numTaxa-1)-l1.csv")
get_output_filepath(netid::String) = joinpath(getDataDir(), "output", "$(netid).csv")
get_estimated_sim_output_filepath(truenet::HybridNetwork) = joinpath(getDataDir(), "est_data_output", "n$((truenet.numTaxa)-1)r$(truenet.numHybrids).csv")
get_estimated_sim_output_filepath(netid::String) = joinpath(getDataDir(), "est_data_output", "$(netid).csv")
copy_csv_template(output_path::String) = cp(joinpath(getDataDir(), "output", "fields.csv"), output_path)
copy_est_csv_template(output_path::String) = cp(joinpath(getDataDir(), "est_data_output", "fields.csv"), output_path)

# DATA LOADING FUNCTIONS
"""
`netid` is just a number if Level1 is true, i.e. "100"
"""
function loadPerfectData(netid::String, replicatenum::Int64, maxsize::Int64, dmethod::String; all_have_outgroup::Bool = false, outgroup_removed_after_reroot::Bool = false, verify_constraints::Bool = true, level1::Bool = false)
    file_path = level1 ? getNetworkFilepathLevel1(parse(Int64, netid)) : getNetworkFilepath(netid)
    truenet = readMultiTopology(file_path)[replicatenum]
    avg_bl = get_avg_bl(truenet)
    newick = writeTopology(truenet)
    newick = "($(newick[1:(length(newick)-1)]):$(avg_bl),OUTGROUP:1.0);"
    truenet = readTopology(newick)

    constraints = sateIdecomp(majorTree(truenet), maxsize)
    if all_have_outgroup
        for (i, c) in enumerate(constraints)
            constraints[i] = union(c, ["OUTGROUP"])
        end
    end
    constraints = pruneTruthFromDecomp(truenet, constraints)
    if outgroup_removed_after_reroot
        for c in constraints
            rootatnode!(c, "OUTGROUP")
            deleteleaf!(c, "OUTGROUP")
        end
    end
    
    # If any constraints have root-retics, remove them
    for c in constraints
        hybridbools = [edge.hybrid for edge in c.node[c.root].edge]
        while any(hybridbools)
            PhyloNetworks.deletehybridedge!(c, c.node[c.root].edge[hybridbools][1])
            hybridbools = [edge.hybrid for edge in c.node[c.root].edge]
        end
    end

    if verify_constraints
        InPhyNet.check_constraints(constraints)
    end

    D, namelist = (nothing, nothing)
    if dmethod == "internode_count" || dmethod == "AGIC"
        D, namelist = majorinternodecount(truenet)
    elseif dmethod == "internode_distance" || dmethod == "AGID"
        D, namelist = majorinternodedistance(truenet)
    else
        error("Unrecognized distance method specified.")
    end
    return truenet, constraints, D, namelist
end


@inline function get_avg_bl(net::HybridNetwork)
    bl_sum = 0.
    num_edges = 0
    for e in net.edge
        if e.length != -1. && e.length != 0.
            bl_sum += e.length
            num_edges += 1
        end
    end
    return bl_sum / num_edges
end


get_all_net_ids() = ["n50r2", "n50r5", "n100r5", "n100r10", "n200r10", "n200r20", "n500r25", "n500r50", "n1000r50", "n1000r100"]


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
    robD3 = robustGauss(truenet, constraints, μ=3., σ=3., nsim=nsim)
    robD4 = robustGauss(truenet, constraints, μ=4., σ=4., nsim=nsim)

    # 3. NNI robustness tests
    a0, b0, c0 = robustNNI(truenet, constraints, repeat([0.5], length(constraints)), nsim=2*nsim)
    a1, b1, c1 = robustNNI(truenet, constraints, repeat([1], length(constraints)), nsim=nsim)
    a2, b2, c2 = robustNNI(truenet, constraints, repeat([2], length(constraints)), nsim=nsim)
    a3, b3, c3 = robustNNI(truenet, constraints, repeat([3], length(constraints)), nsim=nsim)
    a4, b4, c4 = robustNNI(truenet, constraints, repeat([4], length(constraints)), nsim=nsim)
    b0 = sum(b0, dims=1)[1,:]
    b1 = sum(b1, dims=1)[1,:]
    b2 = sum(b2, dims=1)[1,:]
    b3 = sum(b3, dims=1)[1,:]
    b4 = sum(b4, dims=1)[1,:]
    c0 = sum(c0, dims=1)[1,:]
    c1 = sum(c1, dims=1)[1,:]
    c2 = sum(c2, dims=1)[1,:]
    c3 = sum(c3, dims=1)[1,:]
    c4 = sum(c4, dims=1)[1,:]

    r1 = getNetDistances(truenet, mnet)
    r2 = DataFrame(
        gaussMean=[
            repeat([1.], length(robD1));
            repeat([2.], length(robD2));
            repeat([3.], length(robD3));
            repeat([4.], length(robD3))
        ],
        gaussStd=[
            repeat([1.], length(robD1));
            repeat([2.], length(robD2));
            repeat([3.], length(robD3));
            repeat([4.], length(robD3))
        ],
        dists=vcat(robD1, robD2, robD3, robD4)
    )
    r3 = DataFrame(
        nniMoves=[
            repeat([0.5], length(a0));
            repeat([1], length(a1));
            repeat([2], length(a2));
            repeat([3], length(a3));
            repeat([4], length(a4))
        ],
        edgeheights=vcat(c0, c1, c2, c3, c4),
        constraintdists=vcat(b0, b1, b2, b3, b4),
        estdists=vcat(a0, a1, a2, a3, a4)
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