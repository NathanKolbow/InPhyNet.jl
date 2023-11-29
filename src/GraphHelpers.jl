using Graphs, PhyloNetworks

"""
Converts the tree/network `net` into a SimpleGraph to leverage already
implemented pathfinding algorithms.

# Arguments
- includeminoredges (default=true): if true, the entire network is translated to a graph.
      Otherwise, only tree-like edges (other than those in `alwaysinclude`) are retained.
- alwaysinclude (default=nothing): edges that should always be included in the graph,
      regardless of the value of `includeminoredges`
"""
function Graph(net::HybridNetwork; includeminoredges::Bool=true, alwaysinclude::Union{Edge,Nothing}=nothing)
    graph = SimpleGraph(net.numNodes)
    for edge in net.edge
        if !includeminoredges && edge.hybrid && !edge.isMajor && edge != alwaysinclude continue end
        nidx1 = findfirst(net.node .== [edge.node[1]])
        nidx2 = findfirst(net.node .== [edge.node[2]])

        if nidx1 !== nothing && nidx2 !== nothing
            add_edge!(graph, nidx1, nidx2)  # edges are undirected by default, don't need to add twice
        end
    end
    return graph
end


"""

Helper function that removes redundant edges in the graph that exist
when some hybrids are pruned from the graph.
"""
function removeredundantedges!(graph::SimpleGraph; keeproot::Union{Nothing,HybridNetwork}=nothing)
    rootidx = -1
    if keeproot !== nothing
        rootidx = keeproot.root
    end

    redundantidxs = []
    for (i, adjvec) in enumerate(graph.fadjlist)
        if length(adjvec) == 2 && i != rootidx
            push!(redundantidxs, i)
        end
    end

    if length(redundantidxs) == 0 return end
    for idx in redundantidxs
        from = graph.fadjlist[idx][1]
        to = graph.fadjlist[idx][2]
        
        vecreplace!(graph.fadjlist[from], idx, to)
        vecreplace!(graph.fadjlist[to], idx, from)
        graph.fadjlist[idx] = []
    end
    return redundantidxs
end


"""

Helper function that adds edges coming out from the root in `net` back into `graph`.
"""
function addrootedge!(graph::SimpleGraph, net::HybridNetwork)
    totaladded = 0
    rootnode = net.node[net.root]
    rootedges = rootnode.edge

    dst1 = ifelse(rootedges[1].node[1] == rootnode, rootedges[1].node[2], rootedges[1].node[1])
    dst2 = ifelse(rootedges[2].node[1] == rootnode, rootedges[2].node[2], rootedges[2].node[1])
    idxdst1 = findfirst(net.node .== [dst1])
    idxdst2 = findfirst(net.node .== [dst2])
    idxsrc = findfirst(net.node .== [rootnode])

    remove_edge!(graph, idxdst1, idxdst2)
    remove_edge!(graph, idxsrc, idxdst1)
    remove_edge!(graph, idxsrc, idxdst2)

    for edge in net.node[net.root].edge
        nidx1 = findfirst(net.node .== [edge.node[1]])
        nidx2 = findfirst(net.node .== [edge.node[2]])
        
        if !(nidx2 in graph.fadjlist[nidx1])
            totaladded += 1
            add_edge!(graph, nidx1, nidx2)
        end
    end
    return totaladded
end

function remove_edge!(graph::SimpleGraph, idx1::Int64, idx2::Int64)
    remidxs = findall(graph.fadjlist[idx1] .== idx2)
    deleteat!(graph.fadjlist[idx1], remidxs)
    
    remidxs = findall(graph.fadjlist[idx2] .== idx1)
    deleteat!(graph.fadjlist[idx2], remidxs)
end