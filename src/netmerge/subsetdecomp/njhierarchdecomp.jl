using Clustering, DataFrames


function njHierarchDecomp(estgts::Vector{HybridNetwork}, maxtax::Int64; clustkwargs...)
    D, namelist = calculateAGID(estgts)
    njtree = nj(DataFrame(D, namelist))
    return njHierarchDecomp(njtree, maxtax; clustkwargs...)
end


function njHierarchDecomp(tree::HybridNetwork, maxtax::Int64; clustkwargs...)
    D, namelist = internodedistance(tree)
    hres = hclust(D; clustkwargs...)

    nclust = 1
    hclusts = nothing
    maxsize = Inf
    while maxsize > maxtax
        nclust += 1
        hclusts = cutree(hres, k=nclust)
        maxsize = maximum([sum(hclusts .== i) for i=1:nclust])
    end

    subsets = []
    for uqval in unique(hclusts)
        push!(subsets, namelist[hclusts .== uqval])
    end
    if minimum([length(s) for s in subsets]) < 7
        @warn "WARNING: smallest subset has $(minimum([length(s) for s in subsets])) taxa"
    end
    return subsets
end


function pruneTruthFromDecomp(truenet, subsets)
    nets = Vector{HybridNetwork}()
    for set in subsets
        tempnet = readTopology(writeTopology(truenet))
        for leaf in tempnet.leaf
            if !(leaf.name in set)
                PhyloNetworks.deleteleaf!(tempnet, leaf.name)
            end
        end
    end
    return nets
end