# Helper functions for manipulating PhyloNetworks network structures

# TODO: keep track of relevant `ldict` nodes in the IPT structure and then only copy THOSE
#       here, not every single node
function copyldictcontents!(oldnet::HybridNetwork, newnet::HybridNetwork, ldict::LDict)
    for (newnode, oldnode) in zip(newnet.node, oldnet.node)
        if haskey(ldict, oldnode)
            ldict[newnode] = ldict[oldnode]
        end
    end
end

"""
    splitreticulation!(net::HybridNetwork, retic::Node, reticidx::Real, majorline::LineageNode, minorline::LineageNode, ldict::LDict)

Splits `net` into two other networks. The hybrid node `retic` is turned into a non-hybrid edge in each resultant network.
Resultant networks may have tree topologies.

# Arguments
- `net`: the network to be split
- `retic`: the place to split `net` into two
- `reticidx`: the index of `retic`; passed as an argument to save time searching for it in deep copied networks
- `majorline/minorline`: the lineages the will be attributed to the major and minor edges of `retic` after splitting
- `ldict`: lineage dictionary
"""
function splitreticulation!(net::HybridNetwork, retic::Node, reticidx::Real, majorline::LineageNode, minorline::LineageNode, ldict::LDict)
    majoredge = getparentedge(retic)
    minoredge = getparentedgeminor(retic)

    retic.hybrid = false
    retic.gammaz = -1.
    majoredge.gamma = 1.0
    minoredge.gamma = -1.
    majoredge.hybrid = false
    minoredge.hybrid = false
    minoredge.ismajor = true

    retic.name = "s_"*retic.name

    removehybridreference!(net, retic)

    delnode = nothing   # redundant node that crops up when one of the lineages is empty
    # Adjust the topology
    if nlineages(majorline) != 0
        # Make the new node
        newmajornode = retic

        # Detach the minor edge
        detachedge!(newmajornode, minoredge)
        majoredge.isChild1 = majoredge.node[1] == newmajornode
        minoredge.isChild1 = false

        # Update ldict
        ldict[newmajornode] = majorline
    else
        # Set the delnode
        delnode = getparent(majoredge)

        # Delete major edge
        detachedge!(retic, majoredge)
        majoredge.isChild1 = false

        for node in majoredge.node detachedge!(node, majoredge) end
        PhyloNetworks.deleteEdge!(net, majoredge)
    end

    if nlineages(minorline) != 0
        if nlineages(majorline) != 0
            # Make the new node while attaching it to the minor edge
            newminornode = Node(minimum([node.number for node in net.node]) - 1, false, false, [minoredge])

            # Connect minor edge to the new node
            # `retic/newmajornode` has already been disconnected from `minoredge`
            push!(minoredge.node, newminornode)
            minoredge.isChild1 = minoredge.node[1] == newminornode

            # Add minor node to `net`
            PhyloNetworks.pushNode!(net, newminornode)
        else
            # Update edge settings
            minoredge.isChild1 = minoredge.node[1] == retic

            # Major edge removed and handled above, now just need to update ldict
            ldict[retic] = minorline
        end
    else
        # Set the delnode
        delnode = getparent(minoredge)

        # Major edge dealt with, we just need to remove references from minoredge
        detachedge!(retic, minoredge)

        for node in minoredge.node detachedge!(node, minoredge) end
        PhyloNetworks.deleteEdge!(net, minoredge)
    end

    # Remove the redundant node if it exists
    if delnode !== nothing
        length(delnode.edge) == 2 || error("`delnode` has "*string(length(delnode.edge))*" edges.")

        pedge = getparentedge(delnode)
        cedge = getchildedge(delnode)
        cnode = getchild(delnode)
        pnode = getparent(delnode)

        # Keep `pedge`, delete `cedge`
        pedge.length += cedge.length
        delnode.edge = []
        detachedge!(delnode, pedge)
        detachedge!(delnode, cedge)
        detachedge!(cnode, cedge)

        # Remove references from deleted items
        PhyloNetworks.deleteEdge!(net, cedge)
        PhyloNetworks.deleteNode!(net, delnode)

        # Attach `pedge` to `cnode`
        push!(pedge.node, cnode)
        push!(cnode.edge, pedge)

        # Update attributes
        pedge.isChild1 = cnode == pedge.node[1]
    end

    return nothing
end

"""
    detachedge!(node::Node, edge::Edge)

Used after `splitreticulation!`. Before the function, `node` has 2 parent edges
but we want it to only have 1 now; namely not `edge`.
"""
function detachedge!(node::Node, edge::Edge)
    # Remove the node from the edge's node list
    for (i, tempnode) in enumerate(edge.node)
        if tempnode == node
            deleteat!(edge.node, i)
            break
        end
    end

    # Remove the edge from the node's edge list
    for (i, tempedge) in enumerate(node.edge)
        if tempedge == edge
            deleteat!(node.edge, i)
            break
        end
    end
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

function removehybridreference!(net::HybridNetwork, retic::Node)
    filter!(s -> s ≠ retic, net.hybrid)
    net.numhybrids -= 1
end