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
function splitreticulation!(net::HybridNetwork, retic::PhyloNetworks.Node, leftline::LineageNode, rightline::LineageNode, ldict::LDict)
    println("\nNet address: "*repr(UInt64(pointer_from_objref(net))))
    println("retic address: "*repr(UInt64(pointer_from_objref(retic))))
    println()
    leftedge = getparentedge(retic)
    rightedge = getparentedgeminor(retic)

    retic.hybrid = false
    retic.gammaz = -1.

    # TOOD: WARNING: MAY NOT BE UPDATING ALL NECESSARY 
    # ATTRIBUTES TO MAINTAIN TOPOLOGICAL CONSISTENCY
    retic.name = "SPLIT-"*retic.name
    newleftnode = deepcopy(retic)
    newrightnode = deepcopy(retic)

    replacenodeinedge!(leftedge, retic, newleftnode)
    replacenodeinedge!(rightedge, retic, newrightnode)

    removehybridreference!(net, retic)

    if nlineages(leftline) == 0
        # remove the redundant node
        # -- fuse edges
        removenodeandedge!(leftedge, newleftnode)
    elseif nlineages(rightline) == 0
        # remove the redundant node
        # -- fuse edges
        removenodeandedge!(rightedge, newrightnode)
    end

    ldict[newleftnode] = leftline
    ldict[newrightnode] = rightline

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