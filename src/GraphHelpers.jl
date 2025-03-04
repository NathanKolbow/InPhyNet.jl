using Graphs, PhyloNetworks

"""
Converts the tree/network `net` into a SimpleGraph to leverage already
implemented pathfinding algorithms.

# Arguments
- includeminoredges (default=true): if true, the entire network is translated to a graph.
      Otherwise, only tree-like edges (other than those in `alwaysinclude`) are retained.
- alwaysinclude (default=nothing): edges that should always be included in the graph,
      regardless of the value of `includeminoredges`
- withweights (default=false): return a set of weights corresponding to branch lengths as well
- removeredundantedgecost (default=false): if true, finds all the nodes that would be removed via
      `removeredundantedges!` and sets their weights in `W` to 0
"""
function Graph(net::HybridNetwork; includeminoredges::Bool=true, alwaysinclude::Union{Edge,Nothing}=nothing, withweights::Bool=false, minoredgeweight::Float64=1., removeredundantedgecost::Bool=false)
    graph = SimpleGraph(net.numnodes)
    weights = Matrix{Float64}(undef, net.numnodes, net.numnodes)
    weights .= Inf
    nodemap = Dict{Node, Int64}(node => idx for (idx, node) in enumerate(net.node))
    for edge in net.edge
        if !includeminoredges && edge.hybrid && !edge.ismajor && edge != alwaysinclude continue end
        enode1 = edge.node[1]
        enode2 = edge.node[2]
        if haskey(nodemap, enode1) && haskey(nodemap, enode2)
            add_edge!(graph, nodemap[enode1], nodemap[enode2])
            if withweights
                weight = 1
                if edge.hybrid && !edge.ismajor
                    weight = minoredgeweight
                elseif edge.length == -1.
                    weight = 1
                end
                
                weights[nodemap[enode1], nodemap[enode2]] =
                    weights[nodemap[enode2], nodemap[enode1]] = weight
            end
        end
    end

    if removeredundantedgecost && withweights
        redundant_idxs = removeredundantedges!(deepcopy(graph), net)
        for idx in redundant_idxs
            non_inf_idxs = weights[idx, :] .!= Inf

            # Instead of setting to 0, reduce one of the idxs by 1. This way
            # we don't accidentally reduce the cost of 2 edges to 0, but we
            # get the desired result still.
            weights[idx, [findfirst(non_inf_idxs)]] = weights[[findfirst(non_inf_idxs)], idx] .-= 1.
        end
    end

    if withweights
        return graph, weights
    end
    return graph
end


"""

Helper function that finds a path in `net` from `nodei` to `nodej`. If that path has
2 or more reticulations in it, the cost to include a hybrid in the path is incrementally
increased until only 1 reticulation is in the path.
"""
function find_valid_node_path(net::HybridNetwork, nodei::Node, nodej::Node)
    retics_in_path = 3
    minor_edge_weight = 1.01
    graph = W = nodesinpath = edgesinpath = nothing

    while retics_in_path > 1
        minor_edge_weight += 0.50

        # find shortest path from `nodei` to `nodej`
        graph, W = Graph(net, includeminoredges=true, withweights=true, minoredgeweight=minor_edge_weight, removeredundantedgecost=true)

        idxnodei = findfirst(net.node .== [nodei])
        idxnodej = findfirst(net.node .== [nodej])
        edgepath = a_star(graph, idxnodei, idxnodej, W)

        nodesinpath = Array{Node}(undef, length(edgepath)+1)
        edgesinpath = Array{Edge}(undef, length(edgepath))
        for (i, gedge) in enumerate(edgepath)
            srcnode = net.node[gedge.src]
            dstnode = net.node[gedge.dst]

            if i == 1 nodesinpath[1] = net.node[gedge.src] end
            nodesinpath[i+1] = dstnode

            netedge = filter(e -> (srcnode in e.node) && (dstnode in e.node), dstnode.edge)[1]
            edgesinpath[i] = netedge
        end

        retics_in_path = sum([e.hybrid && !e.ismajor for e in edgesinpath])
    end

    return graph, W, nodesinpath, edgesinpath
end


function find_valid_node_path(net::HybridNetwork, namei::String, namej::String)
    nodei = net.node[findfirst(node.name == namei for node in net.node)]
    nodej = net.node[findfirst(node.name == namej for node in net.node)]
    return find_valid_node_path(net, nodei, nodej)
end


function tester_path(net, namei, namej)
    nodei = net.node[findfirst(node.name == namei for node in net.node)]
    nodej = net.node[findfirst(node.name == namej for node in net.node)]
    return find_valid_node_path(net, nodei, nodej)
end


"""

Helper function that removes redundant edges in the graph that exist
when some hybrids are pruned from the graph.
"""
function removeredundantedges!(graph::SimpleGraph, net::HybridNetwork; keeproot::Bool=true, W::Union{Nothing,Matrix{Float64}}=nothing)
    rootidx = -1
    if keeproot
        rootidx = net.rooti
    end

    redundantidxs = []
    for (i, adjvec) in enumerate(graph.fadjlist)
        if length(adjvec) == 2 && i != rootidx
            push!(redundantidxs, i)
        end
    end

    if length(redundantidxs) != 0
        for idx in redundantidxs
            from = graph.fadjlist[idx][1]
            to = graph.fadjlist[idx][2]
            
            vecreplace!(graph.fadjlist[from], idx, to)
            vecreplace!(graph.fadjlist[to], idx, from)

            if W !== nothing
                W[from, to] = W[to, from] = mean([W[from, idx], W[to, idx]])
                W[from, idx] = W[idx, from] = W[to, idx] = W[idx, to] = 0.
            end

            graph.fadjlist[idx] = []
        end
    end

    # Remove any nodes in the graph that (1) only have 1 connected node and (2) are not leaves
    for (vert_idx, adj_list) in enumerate(graph.fadjlist)
        if keeproot && vert_idx == net.rooti continue
        elseif length(adj_list) == 1 && !net.node[vert_idx].leaf
            graph.fadjlist[vert_idx] = []
            graph.fadjlist[adj_list[1]] = setdiff(graph.fadjlist[adj_list[1]], vert_idx)
            push!(redundantidxs, vert_idx)
        end
    end
    
    if length(redundantidxs) == 0
        return []
    else
        # If something was removed, recurse.
        #
        # We do this b/c this method removes redundant edges in 2 ways and either one
        # could lead to the other needing to be run again.
        return vcat(redundantidxs, removeredundantedges!(graph, net, keeproot=keeproot, W=W))
    end
    return redundantidxs
end


"""

Helper function that adds edges coming out from the root in `net` back into `graph`.
"""
function addrootedge!(graph::SimpleGraph, net::HybridNetwork)
    totaladded = 0
    rootnode = getroot(net)
    rootedges = rootnode.edge

    dst1 = ifelse(rootedges[1].node[1] == rootnode, rootedges[1].node[2], rootedges[1].node[1])
    dst2 = ifelse(rootedges[2].node[1] == rootnode, rootedges[2].node[2], rootedges[2].node[1])
    idxdst1 = findfirst(net.node .== [dst1])
    idxdst2 = findfirst(net.node .== [dst2])
    idxsrc = findfirst(net.node .== [rootnode])

    remove_edge!(graph, idxdst1, idxdst2)
    remove_edge!(graph, idxsrc, idxdst1)
    remove_edge!(graph, idxsrc, idxdst2)

    for edge in getroot(net).edge
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