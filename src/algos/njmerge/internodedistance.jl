
function majorinternodedistance(N::HybridNetwork)
    return internodedistance(majorTree(N))
end

function internodedistance(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}}=nothing)
    D = zeros(N.numTaxa, N.numTaxa)
    if namelist == nothing
        namelist = [l.name for l in N.leaf]
    end

    Ngraph = Graph(N)
    removeredundantedges!(Ngraph)
    nodelistidx = [findfirst([n.name == name for n in N.node]) for name in namelist]

    for i=1:(N.numTaxa-1)
        nodenumi = nodelistidx[i]
        nodei = N.node[nodenumi]
        for j=(i+1):N.numTaxa
            nodenumj = nodelistidx[j]
            nodej = N.node[nodenumj]

            D[i, j] = D[j, i] = length(a_star(Ngraph, nodenumi, nodenumj)) - 1
        end
    end
    return D, namelist
end

function calculateAGID(Ns::AbstractVector{HybridNetwork})
    D, namelist = internodedistance(Ns[1])
    for j=2:length(Ns)
        D .+= internodedistance(Ns[j], namelist=namelist)[1]
    end
    return D ./ length(Ns), namelist
end

# TODO: remove these functions from here and `netnjmerge.jl` and place it in an appropriate file
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