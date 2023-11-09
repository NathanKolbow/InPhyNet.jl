struct SubNet
    # TODO: Sketch out how merging subnets is going to work
    #       *in detail* before finalizing and using this struct
    nodes::Vector{Node}
    edges::Vector{Edge}
    # some field for where to add next piece?
    # maybe the last node in `nodes` can just be used as reference?
end

# Fundamental methods; the reason this struct exists
"""

Merges two unrooted SubNets into a single unrooted SubNet
"""
function mergesubnets!(n1::SubNet, n2::SubNet)
    error("Not implemented yet")
    
end


"""

Converts the final SubNet in the net nj merge algo into a
    HybridNetwork object.
"""
function HybridNetwork(subnet::SubNet)
    error("Not implemented")
end