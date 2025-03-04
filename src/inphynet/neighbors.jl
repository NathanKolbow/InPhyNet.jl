using PhyloNetworks


"""
Gets the set of neighbors in network `net`.
Only implemented for trees right now.
"""
function get_neighbor_set(net::HybridNetwork)

    net.numhybrids == 0 || throw(ArgumentError("Only implemented for trees right now."))
    pairs = Set{Tuple{Node, Node}}()

    # Inefficient implementation, but this algo will never be a bottleneck so who cares!
    for leaf in net.leaf
        neighbor_nodes = getnodes(getparent(leaf))
        for node in neighbor_nodes
            if node.leaf && node != leaf
                # This way we don't enter the same pair twice
                if leaf.name < node.name
                    push!(pairs, (leaf, node))
                else
                    push!(pairs, (node, leaf))
                end
            end
        end
    end

    return pairs

end