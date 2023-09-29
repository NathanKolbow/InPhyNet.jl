# Helper functions for manipulating PhyloNetworks network structures

function copyldictcontents!(oldnet::HybridNetwork, newnet::HybridNetwork, ldict::LDict)
    for (newnode, oldnode) in zip(newnet.node, oldnet.node)
        ldict[newnode] = ldict[oldnode]
    end
end

# Splits the given reticulation into 2 tree-like edges
# Creates a deep copy of the net in the process
# `retic` is the hybrid node that 2 edges point into
#
# "left" and "right" don't actually hold any meaning here,
# they are just a useful naming convention to keep things clear
function splitreticulation!(net::HybridNetwork, retic::PhyloNetworks.Node, majorline::LineageNode, minorline::LineageNode, ldict::LDict)
    println("\nNet address: "*repr(UInt64(pointer_from_objref(net))))
    println("retic address: "*repr(UInt64(pointer_from_objref(retic))))
    println()
    majoredge = getparentedge(retic)
    minoredge = getparentedgeminor(retic)

    retic.hybrid = false
    retic.gammaz = -1.
    majoredge.gamma = 1.0
    minoredge.gamma = 1.0

    # TOOD: WARNING: MAY NOT BE UPDATING ALL NECESSARY 
    # ATTRIBUTES TO MAINTAIN TOPOLOGICAL CONSISTENCY
    retic.name = "s_"*retic.name
    newleftnode = deepcopy(retic)
    newrightnode = deepcopy(retic)

    replacenodeinedge!(majoredge, retic, newleftnode)
    replacenodeinedge!(minoredge, retic, newrightnode)

    removehybridreference!(net, retic)

    if nlineages(majorline) == 0
        # remove the redundant node
        # -- fuse edges
        removenodeandedge!(majoredge, newleftnode)
    elseif nlineages(minorline) == 0
        # remove the redundant node
        # -- fuse edges
        removenodeandedge!(minoredge, newrightnode)
    end

    ldict[newleftnode] = majorline
    ldict[newrightnode] = minorline

    return newleftnode, newrightnode
end


function removenodeandedge!(edge::Edge, node::Node)
    # 1. remove `edge` from the other node that it is attached to
    othernode = edge.node[findfirst(edge.node .!== [node])]
    deleteat!(othernode.edge, findfirst(othernode.edge .== [edge]))

    # 2. remove everything from `edge`
    edge.node = []

    # 3. remove everything from `node`
    node.edge = []
end

function replacenodeinedge!(edge::Edge, oldnode::Node, newnode::Node)
    edge.node[findfirst(edge.node .== [oldnode])] = newnode
end

removehybridreference!(net::HybridNetwork, retic::Node) = filter!(s -> s ≠ retic, net.hybrid)