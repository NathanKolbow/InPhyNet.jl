using Distributions, Random, Combinatorics, StatsBase, InPhyNet, PhyloNetworks, StatsBase, DataFrames, CSV, LinearAlgebra


# Main driver for manuscript sim 2(i)
function randomPartitionRobustness(truenet::HybridNetwork, constraintsizes::Vector{Int64}, D, namelist; nsim::Int64=1000)
    esterrors = zeros(nsim) .- 1.
    nhyb = zeros(nsim) .- 1.

    Threads.@threads for iter=1:nsim
        subsets = selectRandomSubsets(constraintsizes, namelist)
        constraints = pruneTruthFromDecomp(truenet, subsets)

        try
            nhyb[iter] = sum([c.numHybrids for c in constraints])

            mnet = netnj(D, constraints, namelist, supressunsampledwarning=true)
            itererror = getNetDistances(truenet, mnet)
            esterrors[iter] = itererror
        catch e
            if typeof(e) != ArgumentError
                @show typeof(e)
                throw(e)
            end
        end

    end

    return esterrors, nhyb
end


# Main driver function for manuscript sim 1(v)
function monophyleticSwappingRobustness(truenet, constraints, D, namelist, nswaps; nsim::Int64=1000)
    if typeof(nswaps) <: AbstractVector nswaps = Vector(nswaps)
    else nswaps = [nswaps] end
    constraintnames = [[l.name for l in c.leaf] for c in constraints]
    
    esterrors = zeros(nsim * length(nswaps)) .- 1.

    for (j, currswaps) in enumerate(nswaps)
        Threads.@threads for iter=1:nsim
            iteridx = (j-1)*nsim + iter

            tempnames = deepcopy(constraintnames)
            for swapiter=1:currswaps
                subsetidxs = sample(1:length(tempnames), 2, replace=false)
                swapidx1 = sample(1:length(tempnames[subsetidxs[1]]), 1)
                swapidx2 = sample(1:length(tempnames[subsetidxs[2]]), 1)

                savename = tempnames[subsetidxs[1]][swapidx1]
                tempnames[subsetidxs[1]][swapidx1] = tempnames[subsetidxs[2]][swapidx2]
                tempnames[subsetidxs[2]][swapidx2] = savename
            end
            constraints = pruneTruthFromDecomp(truenet, tempnames)

            try
                mnet = netnj(D, constraints, namelist, supressunsampledwarning=true)
                itererror = getNetDistances(truenet, mnet)
                esterrors[iteridx] = itererror
            catch e
                if typeof(e) != ArgumentError
                    @show typeof(e)
                    throw(e)
                end
            end
        end
    end

    return esterrors
end


# Main driver function for manuscript sims 1(i)-(iv)
function monophyleticRobustness(truenet::HybridNetwork, constraints::Vector{HybridNetwork}, D::Matrix{<:Real}, namelist::Vector{String}; nsim::Int64=1000)
    totalnnigen = Uniform(0, length(constraints) * 4)

    # Recorded values
    esterrors = zeros(nsim) .- 1.
    constraintdiffs = zeros(length(constraints), nsim) .- 1.
    gausserrors = zeros(nsim) .- 1.
    nretics_est = zeros(nsim) .- 1.
    #

    #
    nrows = size(D, 1)
    std0 = upperTriangStd(D)
    #

    fortime = @elapsed Threads.@threads for iter=1:nsim
        # Randomly generate the Gaussian noise parameters
        gaussMean = gaussSd = rand(Uniform(0, 2*std0))
        gausserrors[iter] = gaussSd

        # Randomly generate the number of NNI moves
        totalnnimoves = Int64(round(rand(totalnnigen)))
        nnimoves = sample(1:length(constraints), totalnnimoves, replace=true)
        nnimoves = Vector{Int64}([sum(nnimoves .== i) for i=1:length(constraints)])

        try
            esterrors[iter], constraintdiffs[:,iter], _, nretics_est[iter] =
                runRobustSim(truenet, constraints, D, namelist, gaussMean, gaussSd, nnimoves)
        catch e
            if typeof(e) != ArgumentError
                @show typeof(e)
                throw(e)
            end
        end
    end
    print("Took $(round(fortime, digits=2)) seconds\n")
    return esterrors, gausserrors, constraintdiffs, nretics_est
end


function upperTriangStd(D::Matrix)
    nrows = size(D, 1)
    idxs = Array{CartesianIndex{2}}(undef, Int64(nrows * (nrows - 1) / 2))
    _idx = 1
    for i = 1:(nrows-1)
        for j = (i+1):nrows
            idxs[_idx] = CartesianIndex(i, j)
            _idx += 1
        end
    end
    return std(D[idxs])
end


function robustNNI(truenet::HybridNetwork, constraints::Vector{HybridNetwork},
    nmoves::Vector{<:Real}; nsim::Int64=100, printruntime::Bool=true)

    # Runtime metrics #
    time_copying = 0.
    time_nni = 0.
    time_constraintRF = 0.
    time_mnet = 0.
    time_mnetRF = 0.
    ###################

    dists = Array{Int64}(undef, nsim)
    constraintdiffs = Array{Int64}(undef, length(constraints), nsim)
    edgeheights = zeros(length(constraints), nsim)

    constraintEdgeHeights = [PhyloNetworks.getHeights(c) for c in constraints]

    Threads.@threads for i=1:nsim
        time_copying += @elapsed newconstraints = copyConstraints(constraints)
        
        for (j, (c, newc, moves)) in enumerate(zip(constraints, newconstraints, nmoves))
            if moves < 1 moves = Int64(rand() < moves) end
            time_nni += @elapsed for _=1:moves
                nniedge = doRandomNNI!(newc)
                edgeheights[j, i] += constraintEdgeHeights[j][findfirst(newc.edge .== [nniedge])] 
            end
            if moves > 0 edgeheights[j, i] /= moves end
            time_constraintRF += @elapsed constraintdiffs[j, i] = hardwiredClusterDistance(newc, c, false)
        end

        time_mnet += @elapsed mnet = runGroundTruthPipeline(truenet, newconstraints)
        @elapsed dists[i] = getNetDistances(truenet, mnet)
    end

    # Runtime metrics #
    if printruntime
        println("Time taken: $(round(time_copying + time_nni + time_constraintRF + time_mnet + time_mnetRF, digits=2))s")
        println("\tcopying: $(round(time_copying, digits=2))s")
        println("\tNNI: $(round(time_nni, digits=2))s")
        println("\tconstraint RF: $(round(time_constraintRF, digits=2))s")
        println("\talgo: $(round(time_mnet, digits=2))s")
        println("\tmerged net RF: $(round(time_mnetRF, digits=2))s")
    end
    ###################

    return dists, constraintdiffs, edgeheights
end


function stackRobustNNI(truenet, constraints, nmovevecs, nsimeach)
    dists = Array{Int64}(undef, length(nmovevecs)*nsimeach)
    constraintdiffs = Array{Int64}(undef, length(constraints), length(nmovevecs)*nsimeach)
    edgeheights = Array{Float64}(undef, length(constraints), length(nmovevecs)*nsimeach)
    
    for j=1:length(nmovevecs)
        sidx = (j-1)*nsimeach + 1
        eidx = j*nsimeach

        dists[sidx:eidx], constraintdiffs[:, sidx:eidx], edgeheights[:, sidx:eidx] = 
            robustNNI(truenet, constraints, nmovevecs[j], nsim=nsimeach, printruntime=false)
    end
    return dists, constraintdiffs, edgeheights
end


function robustGauss(truenet, constraints; μ::Float64=0., σ::Float64=1., nsim::Int64=100)
    dists = Array{Int64}(undef, nsim)
    dists .= -1
    rgen = Normal(μ, σ)
    D, namelist = majorinternodedistance(truenet)

    Threads.@threads for i=1:nsim
        Dcopy = deepcopy(D)
        addnoise!(Dcopy, rgen)

        try
            mnet = netnj(Dcopy, constraints, namelist, supressunsampledwarning=true)
            dists[i] = getNetDistances(truenet, mnet)
        catch e
        end
    end
    return dists[dists .!= -1]
end


function robustUniformProportion(truenet, constraints, p; nsim::Int64=100)
    dists = Array{Int64}(undef, nsim)
    dists .= -1

    Threads.@threads for i=1:nsim
        constraintscopy = copyConstraints(constraints)
        D, namelist = majorinternodedistance(truenet)
        addpnoise!(D, p)
        
        mnet = netnj(D, constraintscopy, namelist, supressunsampledwarning=true)
        dists[i] = getNetDistances(truenet, mnet)
    end
    return dists[dists .!= -1]
end


function robustRandD(truenet, constraints; nsim::Int64=100)
    dists = Array{Int64}(undef, nsim)
    dists .= -1

    Threads.@threads for i=1:nsim
        constraintscopy = copyConstraints(constraints)
        D, namelist = majorinternodedistance(truenet)
        randomize!(D)
        
        mnet = netnj(D, constraintscopy, namelist, supressunsampledwarning=true)
        dists[i] = getNetDistances(truenet, mnet)
    end
    return dists[dists .!= -1]
end


@inline function randomize!(D)
    n = size(D)[1]
    for i=1:n
        for j=(i+1):n
            D[i, j] = D[j, i] = rand()
        end
    end
end


@inline function addpnoise!(D, p)
    n = size(D)[1]
    for i=1:n
        for j=(i+1):n
            val = D[i, j]
            D[i, j] = D[j, i] = rand(Uniform((1-p)*val, (1+p)*val))
        end
    end
end


@inline function addnoise!(D, rgen)
    n = size(D)[1]
    for i=1:n
        for j=(i+1):n
            D[i, j] += rand(rgen)
            D[j, i] = D[i, j]
        end
    end
    if minimum(D) < 0 D .-= minimum(D) end
    for i=1:n D[i,i] = 0 end
end


function doRandomNNI!(net; maxattempts::Int64=100)
    j = 0
    e = sample(net.edge, 1, replace=false)[1]
    while j < 100 && nni!(net, e) === nothing
        e = sample(net.edge, 1, replace=false)[1]
    end
    if j < 100
        return e
    else
        return nothing
    end
end


function getNetDistances(truenet, estnet)
    if "OUTGROUP" in [leaf.name for leaf in truenet.leaf]
        try
            rootatnode!(estnet, "OUTGROUP")
            return hardwiredClusterDistance(truenet, estnet, true)
        catch e
            return hardwiredClusterDistance(truenet, estnet, true)
        end
    elseif truenet.numTaxa > 50
        return hardwiredClusterDistance(truenet, estnet, true)
    else
        return hardwiredClusterDistance(truenet, estnet, false)
    end
end


function runAndSaveRobustnessPipeline(netid::String, whichConstraints::Int64=1)
    truenet, constraints = loadTrueData(netid, whichConstraints)
    fileprefix = "$(netid)-$(whichConstraints)"

    baselineDist = robustDdf = robustNNIdf = nothing
    if !isfile("data/$(fileprefix)_robustDdf.csv")
        totalt = @elapsed baselineDist, robustDdf, robustNNIdf = runGroundTruthRobustnessPipeline(truenet, constraints, nsim=1000)
        println("\nFinished in $(round(totalt/60, digits=2)) minutes.")

        # Save results
        CSV.write("data/$(fileprefix)_robustDdf.csv", robustDdf)
        CSV.write("data/$(fileprefix)_robustNNIdf.csv", robustNNIdf)
        open("data/$(fileprefix)_baselineDist.dat", "w+") do f
            write(f, "$(baselineDist)")
        end
    else
        baselineDist = parse(Float64, readlines("data/$(fileprefix)_baselineDist.dat")[1])
        robustDdf = CSV.read("data/$(fileprefix)_robustDdf.csv", DataFrame)
        robustNNIdf = CSV.read("data/$(fileprefix)_robustNNIdf.csv", DataFrame)
    end
    return baselineDist, robustDdf, robustNNIdf
end


function selectRandomSubsets(constraintsizes::Vector{Int64}, namelist)
    namelist = deepcopy(namelist)
    subsets = Vector{Vector{String}}()

    # Randomly select subsets
    for i=1:length(constraintsizes)
        samp = sample(namelist, constraintsizes[i], replace=false)
        setdiff!(namelist, samp)
        push!(subsets, samp)
    end

    return subsets
end


function runRobustSim(truenet::HybridNetwork, constraints::Vector{HybridNetwork}, D::Matrix{<:Real}, namelist::Vector{String},
    gaussMean::Real, gaussSd::Real, NNImoves::Vector{<:Real}; copyD::Bool=true, copyconstraints::Bool=true)

    if copyD D = deepcopy(D) end
    origconstraints = constraints
    if copyconstraints constraints = copyConstraints(constraints) end

    length(NNImoves) == length(constraints) || error("Must specify same number of NNI moves as exist constraints.")

    # Add noise
    if gaussSd > 0
        addnoise!(D, Normal(gaussMean, gaussSd))
    end

    # Do NNI moves
    constraintdiffs = zeros(length(constraints))
    if maximum(NNImoves) > 0
        for (i, (c, nmoves)) in enumerate(zip(constraints, NNImoves))
            for _=1:nmoves doRandomNNI!(c) end
            constraintdiffs[i] = hardwiredClusterDistance(origconstraints[i], constraints[i], false)
        end
    end
    tempcs = copyConstraints(constraints)

    # Merge the nets
    try
        mnet = netnj(D, constraints, namelist, supressunsampledwarning=true)
        esterror = getNetDistances(truenet, mnet)
        return esterror, constraintdiffs, writeTopology(mnet), mnet.numHybrids
    catch e
        if typeof(e) != ArgumentError
            println("ERROR RECEIVED")
            @show typeof(e)
            println("CONSTRAINTS AFTER NNI MOVES BUT BEFORE MERGING:")
            for c in tempcs println("\t$(writeTopology(c))") end
            println("TRUE NET")
            println("\t$(writeTopology(truenet))")

            throw(e)
        else
            return -1, constraintdiffs, "", -1
        end
    end
end


copyConstraints(cs::Vector{HybridNetwork}) = [readTopology(writeTopology(c)) for c in cs]