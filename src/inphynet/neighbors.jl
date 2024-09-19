using PhyloNetworks


"""
Gets the set of neighbors in network `net`.
Only implemented for trees right now.
"""
function get_neighbor_set(net::HybridNetwork)

    net.numHybrids == 0 || throw(ArgumentError("Only implemented for trees right now."))
    pairs = Set{Tuple{Node, Node}}()

    # Inefficient implementation, but this algo will never be a bottleneck so who cares!
    for leaf in net.leaf
        neighbor_nodes = getnodes(leaf)
        for node in neighbor_nodes
            if node.leaf
                # min/max so that we don't enter the same pair twice
                push!(pairs, (min(leaf, node), max(leaf, node)))
            end
        end
    end

    return pairs

end