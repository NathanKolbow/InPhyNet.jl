struct SubNet
    # TODO: Sketch out how merging subnets is going to work
    #       *in detail* before finalizing and using this struct
    nodes::Vector{Node}
    edges::Vector{Edge}
    # some field for where to add next piece?
    # maybe the last node in `nodes` can just be used as reference?
    SubNet() = error("not yet implemented")
    function SubNet(i, name)
        node = Node(i, true)
        node.name = name
    end
end


"""

Merges two unrooted SubNets into a single unrooted SubNet
IMPORTANT: this function should NOT make copies of edges or
           nodes, object references are stored and used in
           the main merging algorithm
"""
function mergesubnets!(n1::SubNet, n2::SubNet)
    error("Not implemented yet")
    
    return mergednet
end


"""

Converts the final SubNet in the net nj merge algo into a
    HybridNetwork object.
"""
function HybridNetwork(subnet::SubNet)
    error("Not implemented")

    # Remember to re-number edges
end