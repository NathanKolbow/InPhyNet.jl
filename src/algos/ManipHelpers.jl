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
function splitreticulation!(net::HybridNetwork, retic::Node, reticidx::Real, majorline::LineageNode, minorline::LineageNode, ldict::LDict)
    println("\nSPLITTING\nSPLITTING\nSPLITTING\nSPLITTING\nSPLITTING")
    
    majoredge = getparentedge(retic)
    minoredge = getparentedgeminor(retic)

    retic.hybrid = false
    retic.gammaz = -1.
    majoredge.gamma = 1.0
    minoredge.gamma = 1.0

    # TODO: WARNING: MAY NOT BE UPDATING ALL NECESSARY 
    # ATTRIBUTES TO MAINTAIN TOPOLOGICAL CONSISTENCY
    retic.name = "s_"*retic.name
    newleftnode = retic

    removehybridreference!(net, retic)
    if nlineages(majorline) != 0
        println("majorline processed")

        # Swap this node for retic on the major edge
        replacenodeinedge!(majoredge, retic, newleftnode)
        net.node[reticidx] = newleftnode
        
        # Add to lineage dict
        ldict[newleftnode] = majorline
    end
    if nlineages(minorline) != 0
        println("minorline processed")
        
        # Copy the retic node
        newrightnode = deepcopy(retic)

        # Swap this node for retic on the minor edge
        replacenodeinedge!(minoredge, retic, newrightnode)

        # If majorline is empty, some housekeeping
        if nlineages(majorline) == 0
            net.node[reticidx] = newrightnode
        else
            newrightnode.number = minimum([node.number for node in net.node]) - 1
            push!(net.node, newrightnode)
        end
        
        # Add to lineage dict
        ldict[newrightnode] = minorline
    end



    # if nlineages(majorline) == 0
    #     # remove the redundant node
    #     # -- fuse edges
    #     newleftnode.name = "REMOVED"*newleftnode.name
    #     removenodeandedge!(majoredge, newleftnode)
    #     deleteat!(net.node, findfirst([net.node .== [newleftnode]]))
    # elseif nlineages(minorline) == 0
    #     # remove the redundant node
    #     # -- fuse edges
    #     newrightnode.name = "REMOVED"*newrightnode.name
    #     removenodeandedge!(minoredge, newrightnode)
    #     deleteat!(net.node, findfirst([net.node .== [newrightnode]]))
    # end

    # println("\nNet address: "*repr(UInt64(pointer_from_objref(net))))
    # println("retic address: "*repr(UInt64(pointer_from_objref(retic))))
    # println("newleftnode address: "*repr(UInt64(pointer_from_objref(newleftnode))))
    # println("newrightnode address: "*repr(UInt64(pointer_from_objref(newrightnode))))
    # println()

    return nothing
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