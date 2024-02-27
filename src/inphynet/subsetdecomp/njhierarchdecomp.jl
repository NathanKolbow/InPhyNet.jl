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
    nets = Array{HybridNetwork}(undef, length(subsets))
    for (i, set) in enumerate(subsets)
        tempnet = deepcopy(truenet)
        namelist = [leaf.name for leaf in tempnet.leaf]
        for name in namelist
            if !(name in set)
                deleteleaf!(tempnet, name)
            end
        end

        # Sometimes, the given leafset has incoming reticulations that originate
        # from outside of its subnetwork. In this case, the reticulations are kept
        # in the pruned network but given the root as their origin, which is not desired.
        true_dict = Dict(n.number => n for n in truenet.node)
        for temp_node in tempnet.node
            true_node = true_dict[temp_node.number]
            temp_numbers = [n.number for n in getchildren(temp_node)]
            true_numbers = [n.number for n in getchildren(true_node)]

            length(temp_numbers) == length(true_numbers) || continue

            if length(temp_numbers) == 2
                if !((temp_numbers[1] == true_numbers[1] && temp_numbers[2] == true_numbers[2]) ||
                    (temp_numbers[1] == true_numbers[2] && temp_numbers[2] == true_numbers[1]))
                    # Remove the hybrid edges attached here
                    @warn "`pruneTruthFromDecomp` HAS NOT BEEN SUFFICIENTLY TESTED TO TRUST IN ITS OUTPUT YET. VERIFY SOME OF THE RESULTS RECEIVED BEFORE CONTINUING."
                    for temp_edge in temp_node.edge
                        if temp_edge.hybrid PhyloNetworks.deletehybridedge!(tempnet, temp_edge) end
                    end
                end
            elseif length(temp_numbers) == 1
                if temp_numbers[1] != true_numbers[1]
                    # Remove the hybrid edges attached here
                    @warn "`pruneTruthFromDecomp` HAS NOT BEEN SUFFICIENTLY TESTED TO TRUST IN ITS OUTPUT YET. VERIFY SOME OF THE RESULTS RECEIVED BEFORE CONTINUING."
                    for temp_edge in temp_node.edge
                        if temp_edge.hybrid PhyloNetworks.deletehybridedge!(tempnet, temp_edge) end
                    end
                end
            end
        end

        nets[i] = tempnet
    end
    return nets
end