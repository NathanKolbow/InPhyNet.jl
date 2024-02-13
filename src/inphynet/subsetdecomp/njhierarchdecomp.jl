using Clustering, DataFrames


function njHierarchDecomp(estgts::Vector{HybridNetwork}, maxtax::Int64; clustkwargs...)
    D, namelist = calculateAGID(estgts)
    njtree = nj(DataFrame(D, namelist))
    return njHierarchDecomp(njtree, maxtax; clustkwargs...)
end


# TODO: after gathering the clusters, for all pairs of clusters (clusti, clustj)
#       with length(clusti) + length(clustj) < maxtax,
#       select min_(i,j) [ avg_pairwise_dist(clusti, clustj)] using `D`
function njHierarchDecomp(tree::HybridNetwork, maxtax::Int64; clustkwargs...)
    @warn "Function `njHierarchDecomp` is deprecated. Use `sateIdecomp` instead."

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


"""
Returns a list of networks where each entry is a copy of `truenet` pruned to include
only the taxa named in each respective set in `subsets`.

# Example

```julia
mynetwork = readTopology("((t1,t2),(t3,t4));")
pruneTruthFromDecomp(mynetwork, [["t1", "t2"], ["t3", "t4"]])
> 2-element Vector{HybridNetwork}: (t1, t2); and (t3, t4);
```
"""
function pruneTruthFromDecomp(truenet::HybridNetwork, subsets::AbstractVector{<:AbstractVector{<:AbstractString}})
    nets = Vector{HybridNetwork}()
    for set in subsets
        tempnet = readTopology(writeTopology(truenet))
        namelist = [leaf.name for leaf in tempnet.leaf]
        for name in namelist
            if !(name in set)
                PhyloNetworks.deleteleaf!(tempnet, name)
            end
        end
        push!(nets, tempnet)
    end
    return nets
end