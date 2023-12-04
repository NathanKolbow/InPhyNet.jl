using Distributions, Random


function robustNNI(truenet::HybridNetwork, constraints::Vector{HybridNetwork},
    nmoves::Int64; nsim::Int64=100)

    dists = Array{Int64}(undef, nsim)
    constraintdiffs = Array{Int64}(undef, 4, nsim)

    Threads.@threads for i=1:nsim
        newconstraints = copyConstraints(constraints)
        for (j, (c, newc)) in enumerate(zip(constraints, newconstraints))
            for _=1:nmoves
                e = sample(newc.edge, 1, replace=false)[1]
                nni!(newc, e)
            end
            constraintdiffs[j, i] = hardwiredClusterDistance(newc, c, false)
        end

        mnet = runGroundTruthPipeline(truenet, newconstraints)
        dists[i] = hardwiredClusterDistance(truenet, mnet, false)
    end
    return dists, constraintdiffs
end


function robustGauss(truenet, constraints; μ::Float64=0., σ::Float64=1., nsim::Int64=100)
    dists = Array{Int64}(undef, nsim)
    dists .= -1
    rgen = Normal(μ, σ)
    D, namelist = majorinternodedistance(truenet)

    Threads.@threads for i=1:nsim
        constraintscopy = copyConstraints(constraints)
        Dcopy = deepcopy(D)
        addnoise!(Dcopy, rgen)
        
        try
            mnet = netnj(Dcopy, constraintscopy, namelist)
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


copyConstraints(cs::Vector{HybridNetwork}) = [readTopology(writeTopology(c)) for c in cs]