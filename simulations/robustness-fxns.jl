using Distributions, Random


function robustNNI(truenet::HybridNetwork, constraints::Vector{HybridNetwork},
    nmoves::Vector{Int64}; nsim::Int64=100, printruntime::Bool=true)

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

    for i=1:nsim
        time_copying += @elapsed newconstraints = copyConstraints(constraints)
        
        for (j, (c, newc, moves)) in enumerate(zip(constraints, newconstraints, nmoves))
            time_nni += @elapsed for _=1:moves
                nniedge = doRandomNNI!(newc)
                edgeheights[j, i] += constraintEdgeHeights[j][findfirst(newc.edge .== [nniedge])] 
            end
            if moves > 0 edgeheights[j, i] /= moves end
            time_constraintRF += @elapsed constraintdiffs[j, i] = hardwiredClusterDistance(newc, c, false)
        end

        time_mnet += @elapsed mnet = runGroundTruthPipeline(truenet, newconstraints)
        if truenet.numTaxa > 50
            time_mnetRF += @elapsed dists[i] = hardwiredClusterDistance(truenet, mnet, true)
        else
            time_mnetRF += @elapsed dists[i] = hardwiredClusterDistance(truenet, mnet, false)
        end
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
            mnet = netnj(Dcopy, constraints, namelist)
            dists[i] = hardwiredClusterDistance(truenet, mnet, false)
        catch e
            # pass
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
        
        try
            mnet = netnj(D, constraintscopy, namelist)
            dists[i] = hardwiredClusterDistance(truenet, mnet, false)
        catch e
            # pass
        end
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
        
        try
            mnet = netnj(D, constraintscopy, namelist)
            dists[i] = hardwiredClusterDistance(truenet, mnet, false)
        catch e
            # pass
        end
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
            D[i, j] += max(0, rand(rgen))
            D[j, i] = D[i, j]
        end
    end
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


copyConstraints(cs::Vector{HybridNetwork}) = [readTopology(writeTopology(c)) for c in cs]