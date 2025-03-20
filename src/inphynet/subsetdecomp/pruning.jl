"""
Returns a list of networks where each entry is a copy of `net` pruned to include
only the taxa named in each subset contained in `subsets`. `subsets` can also be
an individual vector of strings instead of a nested vector (see examples below).

# Example

```julia
net = readnewick("((t1,t2),(t3,t4));")
prune_network(net, [["t1", "t2"], ["t3", "t4"]])
> 2-element Vector{HybridNetwork}: (t1, t2); and (t3, t4);

prune_network(net, ["t1", "t2", "t3"])
```
"""
function prune_network(net::HybridNetwork, subsets::AbstractVector{<:AbstractVector{<:AbstractString}})
    nets = Array{HybridNetwork}(undef, length(subsets))
    for (i, set) in enumerate(subsets)
        tempnet = deepcopy(net)
        namelist = [leaf.name for leaf in tempnet.leaf]
        for name in namelist
            if !(name in set)
                deleteleaf!(tempnet, name)
            end
        end

        # Sometimes, the given leafset has incoming reticulations that originate
        # from outside of its subnetwork. In this case, the reticulations are kept
        # in the pruned network but given the root as their origin, which is not desired.
        looping = true
        while looping
            looping = false
            true_dict = Dict(n.number => n for n in net.node)
            for temp_node in tempnet.node
                true_node = true_dict[temp_node.number]
                temp_numbers = [n.number for n in getchildren(temp_node)]
                true_numbers = [n.number for n in getchildren(true_node)]

                length(temp_numbers) == length(true_numbers) || continue

                if length(temp_numbers) == 2
                    if !((temp_numbers[1] == true_numbers[1] && temp_numbers[2] == true_numbers[2]) ||
                        (temp_numbers[1] == true_numbers[2] && temp_numbers[2] == true_numbers[1]))
                        # Remove the hybrid edges attached here
                        # @warn "`prune_network` HAS NOT BEEN SUFFICIENTLY TESTED TO TRUST IN ITS OUTPUT YET. VERIFY SOME OF THE RESULTS RECEIVED BEFORE CONTINUING."
                        for temp_edge in temp_node.edge
                            if temp_edge.hybrid
                                PhyloNetworks.deletehybridedge!(tempnet, temp_edge)
                                looping = true
                            end
                        end
                    end
                elseif length(temp_numbers) == 1
                    if temp_numbers[1] != true_numbers[1]
                        # Remove the hybrid edges attached here
                        # @warn "`prune_network` HAS NOT BEEN SUFFICIENTLY TESTED TO TRUST IN ITS OUTPUT YET. VERIFY SOME OF THE RESULTS RECEIVED BEFORE CONTINUING."
                        for temp_edge in temp_node.edge
                            if temp_edge.hybrid
                                PhyloNetworks.deletehybridedge!(tempnet, temp_edge)
                                looping = true
                            end
                        end
                    end
                end
            end
        end

        nets[i] = readnewick(writenewick(tempnet))
    end
    return nets
end
prune_network(net::HybridNetwork, subset::AbstractVector{<:AbstractString}) = 
    prune_network(net, [subset])[1]
prune_networks(nets::AbstractArray{HybridNetwork}, subset::AbstractVector{<:AbstractString}) =
    [prune_network(net, subset) for net in nets]