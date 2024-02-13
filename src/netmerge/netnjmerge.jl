# Author: Nathan Kolbow
# Hosted at https://github.com/NathanKolbow/network-merging
#
# Source code for main components of the network merging algorithm.

function netnj(estgts::Vector{HybridNetwork})
    error("not implemented yet")
end


"""
    netnj(D::Matrix{Float64}, constraints::Vector{HybridNetwork})

Creates a super-network based on the constraint trees/networks in `constraints` and
distance matrix `D`.

# Arguments
- D: distance matrix relating pairs of taxa. This can be generated from estimated gene trees with [`calculateAGID`](@ref)
"""
function netnj(D::Matrix{Float64}, constraints::Vector{HybridNetwork}, namelist::AbstractVector{<:AbstractString}; kwargs...)
    return netnj!(deepcopy(D), Vector{HybridNetwork}(deepcopy(constraints)), namelist)
end

function netnj!(D::Matrix{Float64}, constraints::Vector{HybridNetwork}, namelist::AbstractVector{<:AbstractString}; supressunsampledwarning=false)
    
    PhyloNetworks.check_distance_matrix(D)
    check_constraints(constraints, false)
    n = size(D, 1)

    # If namelist is not provided
    if isempty(namelist)
        namelist = string.(1:n)
    end
    if length(namelist) != n
        m = length(namelist)
        throw(ArgumentError("D has dimensions $n x $n but only $m names were provided."))
    elseif length(unique(namelist)) != length(namelist)
        throw(ArgumentError("Names must be unique."))
    end

    # Empty network
    subnets = Vector{SubNet}([SubNet(i, namelist[i]) for i in 1:n])
    reticmap = ReticMap(constraints)

    # Edge case: remove (but log) retics that emerge from the root of constraints
    rootretics = Array{Union{Edge, Nothing}}(undef, length(subnets))
    for (i, c) in enumerate(constraints)
        hybridbools = [edge.hybrid for edge in c.node[c.root].edge]
        if any(hybridbools)
            sum(hybridbools) == 1 || error("Two reticulations coming out of root, this is not accounted for!")
            hybedge = c.node[c.root].edge[hybridbools][1]
            rootretics[i] = hybedge

            if !supressunsampledwarning
                @warn "Hybridization involving upsampled taxa detected in constraint network $(i). "*
                    "Merging may not behave as expected, see this post for important details: "*
                    "POST NOT MADE YET - PLEASE POST AN ISSUE ON GITHUB"
            end
        else
            rootretics[i] = nothing
        end
    end
    rootreticprocessed = [false for _ in 1:length(constraints)]

    # Main algorithm loop
    while n > 1
        # DEBUG STATEMENT
        # @show n
        possible_siblings = findvalidpairs(constraints, namelist)
        
        # Find optimal (i, j) idx pair for matrix Q
        i, j = findoptQidx(D, possible_siblings)

        # connect subnets i and j
        subnets[i], edgei, edgej = mergesubnets!(subnets[i], subnets[j])
        updateconstraints!(namelist[i], namelist[j], constraints, reticmap, edgei, edgej)

        # if a constraint with a root-retic is down to a single taxa, place that root-retic in the appropriate subnet
        for (cidx, c) in enumerate(constraints)
            if length(c.leaf) == 1 && rootretics[cidx] !== nothing && !rootreticprocessed[cidx]
                fakesubnet = SubNet(-cidx, "__placeholder_for_rootretic_num$(cidx)__")
                subnets[i], edgei, _ = mergesubnets!(subnets[i], fakesubnet)

                try
                    logretic!(reticmap, rootretics[cidx], edgei, "from")
                catch e
                    logretic!(reticmap, rootretics[cidx], edgei, "to")
                end
                # println("logging retic for edge number $(rootretics[cidx].number)")
                rootreticprocessed[cidx] = true
            end
        end

        # collapse taxa i into j
        for l in 1:n
            if l != i && l != j
                D[l, i] = D[i, l] = (D[l, i] + D[j, l] - D[i, j]) / 2
                
                # try this triangle equality updating rule instead
                # a = D[l, i]
                # b = D[l, j]
                # c = D[i, j]
                # # J = acos(-(a^2-b^2-c^2)/(2*b*c))
                # # d = sqrt(b^2 + (c/2)^2 - b*c*cos(J))
                # d = sqrt(b^2/2 + a^2/2)

                # D[l, i] = D[i, l] = d
            end
        end

        # Remove data elements that corresponded to `j`
        idxfilter = [1:(j-1); (j+1):n]
        D = view(D, idxfilter, idxfilter)   # remove j from D
        subnets = view(subnets, idxfilter)
        namelist = view(namelist, idxfilter)

        n -= 1
    end

    mnet = HybridNetwork(subnets[1].nodes, subnets[1].edges)
    mnet.root = mnet.numNodes
    mnet.node[mnet.root].name = "root"

    mnet = placeretics!(mnet, reticmap)
    removeplaceholdernames!(mnet)

    return mnet
end


@inline function removeplaceholdernames!(mnet::HybridNetwork)
    for node in mnet.node
        if node.leaf && startswith(node.name, "__placeholder_for_rootretic_num")
            deleteleaf!(mnet, node)
        elseif !node.leaf && startswith(node.name, "___internal")
            node.name = ""
        end
    end
end


"""

Checks validity of input constraints. So far, the only check is to make sure
that all nodes have exactly 3 edges except for the root.
"""
function check_constraints(constraints::Vector{HybridNetwork}, requirerooted::Bool)
    for net in constraints
        if net.numTaxa == 1 continue end
        for node in net.node
            nedge = length(node.edge)
            if node.leaf
                if length(node.edge) != 1
                    throw(ArgumentError("Leaf nodes must have exactly 1 attached edge."))
                end
            elseif node == net.node[net.root]
                if length(node.edge) != 2 && requirerooted
                    throw(ArgumentError("Root node must have exactly 2 attached edges."))
                end
            elseif length(node.edge) != 3
                throw(ArgumentError("Internal nodes must have exactly 3 attached edges."))
            end
        end
    end
end


"""

Updates `net` to include the reticulations that we've kept track of along the way
in our algo but haven't placed yet.
"""
function placeretics!(net::HybridNetwork, reticmap::ReticMap)
    namepairs = []
    counter = 0
    
    check_reticmap(reticmap)
    for retic in keys(reticmap.map)
        from = reticmap.map[retic][1]
        to = reticmap.map[retic][2]
        
        if from.node[1].name == ""
            from.node[1].name = "___internal"*string(counter+=1)
        end
        if from.node[2].name == ""
            from.node[2].name = "___internal"*string(counter+=1)
        end
        
        if to.node[1].name == ""
            to.node[1].name = "___internal"*string(counter+=1)
        end
        if to.node[2].name == ""
            to.node[2].name = "___internal"*string(counter+=1)
        end

        # Alphanumeric name ordering for easier uniqueness checking
        name1a, name1b = from.node[1].name < from.node[2].name ? [from.node[1].name, from.node[2].name] : [from.node[2].name, from.node[1].name]
        name2a, name2b = to.node[1].name < to.node[2].name ? [to.node[1].name, to.node[2].name] : [to.node[2].name, to.node[1].name]
        newpair = ([name1a, name1b], [name2a, name2b])

        if !(newpair in namepairs)
            push!(namepairs, newpair)
        end
    end
    mnet = readTopology(writeTopology(net))

    edgepairs = []
    for ((fromname1, fromname2), (toname1, toname2)) in namepairs
        fromedge = nothing
        toedge = nothing

        for edge in mnet.edge
            nodenames = [n.name for n in edge.node]
            if fromname1 in nodenames && fromname2 in nodenames
                fromedge = edge
            elseif toname1 in nodenames && toname2 in nodenames
                toedge = edge
            end
        end
        if fromedge === nothing
            error("from edge nothing")
        elseif toedge === nothing
            error("to edge nothing")
        end
        push!(edgepairs, [fromedge, toedge])
    end

    i = 1
    while i <= length(edgepairs)
        fromedge, toedge = edgepairs[i]

        hybnode, hybedge = addhybridedge!(mnet, fromedge, toedge, true)
        mnet.root = findfirst([n.name == "root" for n in mnet.node])

        for j=(i+1):length(edgepairs)
            if edgepairs[2] == toedge
                edgepairs[2] = hybedge
            end
        end
        i += 1
    end

    return mnet
end


function namepairsreplace!(namepairs::AbstractVector, newname::String, oldname1::String, oldname2::String)
    for (_, topair) in namepairs
        if topair[1] == oldname1 || topair[1] == oldname2
            topair[1] = newname
        end
        if topair[2] == oldname1 || topair[2] == oldname2
            topair[2] = newname
        end
    end
end


"""

Updates constraint networks after (i, j) with names (nodenamei, nodenamej) have been
merged. Also updates `reticmap` to keep track of any reticulations that get removed
in this process.

By convention we keep `nodenamei` and replace node names with `nodenamej`
"""
function updateconstraints!(nodenamei::AbstractString, nodenamej::AbstractString, 
    constraints::Vector{HybridNetwork}, reticmap::ReticMap, subnetedgei::Edge, subnetedgej::Edge)

    for (netidx, net) in enumerate(constraints)
        # DEBUG STATEMENTS
        # println("\t$(netidx): $(writeTopology(net))")
        idxi = -1
        idxj = -1
        for (nodeidx, node) in enumerate(net.node)
            if node.name == nodenamei idxi = nodeidx
            elseif node.name == nodenamej idxj = nodeidx
            end
        end

        hasi = idxi != -1
        hasj = idxj != -1
        # if hasi && !hasj:  do nothing
        # if !hasi && !hasj: do nothing
        if hasj && !hasi
            net.node[idxj].name = nodenamei
        elseif hasi && hasj
            nodei = net.node[idxi]
            nodej = net.node[idxj]

            mergeconstraintnodes!(net, nodei, nodej, reticmap, subnetedgei, subnetedgej)
        end
    end
end


function mergeconstraintnodes!(net::HybridNetwork, nodei::Node, nodej::Node, reticmap::ReticMap, subnetedgei::Edge, subnetedgej::Edge)
    parentsi = getnodes(nodei)
    parentsj = getnodes(nodej)

    length(parentsi) == 1 || length(net.node) == 1 || error("Found $(length(parentsi)) nodes above a leaf?")    # sanity check, remove when finalized
    length(parentsj) == 1 || length(net.node) == 1 || error("Found $(length(parentsj)) nodes above a leaf?")    # sanity check; remove when finalized
    
    parentsi = parentsi[1]
    parentsj = parentsj[1]

    if (parentsi == parentsj && length(net.leaf) == 2) || (parentsi == nodej && parentsj == nodei)
        # println("a: ($(nodei.name), $(nodej.name))")

        # TODO: clean this up, when they're nothing we're assigning them randomly right now
        for edge in net.edge
            if edge.hybrid && !edge.isMajor
                if reticmap.map[edge][1] === nothing
                    logretic!(reticmap, edge, subnetedgei, "from")
                end
                if reticmap.map[edge][2] === nothing
                    logretic!(reticmap, edge, subnetedgej, "to")
                end
            end
        end
        ################

        for edge in net.edge deleteEdge!(net, edge) end
        deleteNode!(net, nodej)
        net.node = [nodei]
        net.edge = []
        net.numTaxa = 1
        net.numNodes = 1
        net.numEdges = 0
        net.root = 1
        nodei.edge = []
    elseif parentsi == parentsj
        # println("b: ($(nodei.name), $(nodej.name))")

        # no reticulations: just merge the nodes
        parent = parentsi
        
        # 1. remove the edge connecting nodej and its parent
        parent.edge = filter(e -> e != nodej.edge[1], parent.edge)
        deleteEdge!(net, nodej.edge[1], part=false)
        nodej.edge = []

        # 2. remove node j entirely
        deleteNode!(net, nodej)

        internal = [node for node in getnodes(parent) if node != nodei && node != nodej][1]
        edge1 = nodei.edge[1]
        edge2 = ifelse(parent.edge[1] == edge1, parent.edge[2], parent.edge[1])

        # 3. the links are now: (internal) --edge2-- (parent) --edge1-- (nodei)
        #              we want: (internal) --edge2-- (nodei)
        #    a) remove parent
        #    b) remove edge1
        #    c) remove edge2
        #    c) connect internal and nodei w/ a new edge
        parent.edge = []
        deleteNode!(net, parent)
        
        deleteEdge!(net, edge1)
        deleteEdge!(net, edge2)
        internal.edge = filter(e -> e != edge1 && e != edge2, internal.edge)
        nodei.edge = filter(e -> e != edge1 && e != edge2, nodei.edge)

        newedge = connectnodes!(nodei, internal)  # handy fxn from SubNet.jl
        push!(net.edge, newedge)
    else
        # println("c: ($(nodei.name), $(nodej.name))")
        # println("Before any operations:\n\t$(writeTopology(net))")

        # find shortest path from `nodei` to `nodej`
        graph = Graph(net, includeminoredges=true)

        idxnodei = findfirst(net.node .== [nodei])
        idxnodej = findfirst(net.node .== [nodej])
        edgepath = a_star(graph, idxnodei, idxnodej)

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

        # find the node that should be the new tip after merging
        newtip = nothing
        for node in nodesinpath
            if node.leaf continue end
            isnewtip = true
            for edge in node.edge
                if edge.hybrid && !edge.isMajor
                    isnewtip = false
                    break
                end
            end
            if isnewtip
                newtip = node
                break
            end
        end
        # @show newtip.name
        # @show [n.name for n in nodesinpath]
        
        ishybedge = [e.hybrid && !e.isMajor for e in edgesinpath]
        hybedge = nothing
        hybedgelogged = false
        if newtip === nothing && sum(ishybedge) == 1
            # we follow a hybrid at some point; the new tip is the shared parent of the node at the beginning
            # and the node at the end of the hybrid edge
            hybedge = edgesinpath[ishybedge][1]

            if getparent(hybedge.node[1]) == getparent(hybedge.node[2])
                newtip = getparent(hybedge.node[1])
                for edge in newtip.edge
                    if hybedge.node[1] in edge.node
                        push!(edgesinpath, edge)
                    elseif hybedge.node[2] in edge.node
                        push!(edgesinpath, edge)
                    end
                end
            else
                if length(getparents(hybedge.node[1])) == 2
                    newtip = hybedge.node[2]
                    for e in hybedge.node[1].edge
                        if !(e in edgesinpath)
                            push!(edgesinpath, e)
                        end
                    end
                else
                    newtip = hybedge.node[1]
                    for e in hybedge.node[2].edge
                        if !(e in edgesinpath)
                            push!(edgesinpath, e)
                        end
                    end
                end
            end

            # direction shouldn't lead to any rooting issues b/c we're at the leaves,
            # but direction not being retained here is bad
            hybedgelogged = true
            logretic!(reticmap, hybedge, subnetedgei, "from")
            logretic!(reticmap, hybedge, subnetedgej, "to")
        elseif newtip === nothing && any(ishybedge)
            error("Unknown case. Multiple hybrid edges in the path.")
        elseif newtip === nothing
            println(net)
            println(nodei)
            println(nodej)
            error("Unknown case. No newtip found but we also don't take a hybrid path.")
        end

        # find the retics that we need to keep track of
        # our a_star path goes from `i` to `j`, so any retics we find are relevant to `i` up
        # until we cross the new tip, at which point they become relevant to `j`
        relevanttoi = true
        for node in nodesinpath
            if node == newtip
                relevanttoi = false
            end
            for edge in node.edge
                if edge.hybrid && !edge.isMajor
                    fromorto = ifelse(getChild(edge) == node, "to", "from")
                    if subnetedgei == subnetedgej
                        @error("equiv edges C")
                    end

                    if edge != hybedge || !hybedgelogged
                        logretic!(reticmap, edge, ifelse(relevanttoi, subnetedgei, subnetedgej), fromorto)
                    end
                end
            end
        end

        # purge all the nodes we've passed through to this point from the network
        for node in nodesinpath
            if node == newtip continue end
            deleteNode!(net, node)
        end
        for edge in edgesinpath
            deleteEdge!(net, edge)
            for node in edge.node
                node.edge = filter(e -> e != edge, node.edge)
            end
        end

        newtip.leaf = true
        push!(net.leaf, newtip)
        net.numTaxa += 1
        newtip.name = nodei.name

        fuseredundantedges!(net)
    end
end


"""

Helper function to remove redundant edges (e.g. A --> (internal node) --> (internal node))
that arise from one case of the function `mergeconstraintnodes!`.

Redundant nodes in this case will always have only two edges while *not* being the root.
"""
function fuseredundantedges!(net::HybridNetwork)
    for (i, node) in enumerate(net.node)
        if node.leaf || node.hybrid || node == net.node[net.root] continue end
        if length(node.edge) == 2
            fuseedgesat!(i, net)
        end
    end
end


"""
    findoptQ(D::Matrix{Float64}, idxpairs::Vector{Tuple{<:Integer, <:Integer}})

Finds the minimizer (i*, j*) among all pairs (i, j) in idxpairs for Q, a matrix computed from D.
"""
TIEWARNING = false
function findoptQidx(D::AbstractMatrix{Float64}, validpairs::Matrix{<:Integer})
    global TIEWARNING

    idxpairs = reduce(vcat, [[(i, j) for j in (i+1):size(D,1) if validpairs[i,j] == 1] for i in 1:size(D, 1)])
    if length(idxpairs) == 0
        throw(ArgumentError("No valid idx pairs received in `findoptQidx`."))
    end

    n = size(D)[1]
    sums = sum(D, dims=1)
    best = Inf
    minidx = (-1, -1)

    for (i, j) in idxpairs
        qij = (n-2) * D[i,j] - sums[i] - sums[j]
        if qij < best
            best = qij
            minidx = (i, j)
        elseif qij == best && !TIEWARNING
            TIEWARNING = true
            @warn "Found a tie in the distance matrix, result may not be unique. This warning is shown once per julia session."
        end
    end
    
    return minidx
end


"""
    findvalidpairs(constraints::Vector{HybridNetwork}, namelist::AbstractVector{<:AbstractString})

Finds all valid sibling pairs among the constraint networks.
"""
function findvalidpairs(constraints::Vector{HybridNetwork}, namelist::AbstractVector{<:AbstractString})
    n = length(namelist)

    # Shorthand functions for name lookups (we'll be doing a lot of these if there are many constraints)
    namedict = Dict{AbstractString, Int64}([name => i for (i, name) in enumerate(namelist)])
    idx(name::AbstractString) = namedict[name]

    # initialize matrix
    validpairs = Matrix{Int64}(undef, n, n)     # -1: pair not seen together yet
                                                #  0: pair invalid
                                                #  1: pair valid
    validpairs .= -1

    # go through the constraint networks and validate/invalidate pairs
    for net in constraints
        if net.numTaxa == 1 continue end
        leafidxs = [idx(leaf.name) for leaf in net.leaf]

        # For each pair, override 1 => -1 so that we can invalidate pairs that are not seen together in this net
        netpairs = view(validpairs, leafidxs, leafidxs)
        netpairs[netpairs .== 1] .= -1

        # Find valid sibling pairs
        nodepairs = findsiblingpairs(net)   # returned as nodes, need to convert to idxs
        nodestoidx(nodepair) = CartesianIndex(idx(nodepair[1].name), idx(nodepair[2].name))
        pairidxs = map(nodestoidx, nodepairs)

        # Set valid pair idxs to 1 only if they are either 1 or -1
        # if a pair idx is 0 then it was invalid elsewhere and needs to stay 0
        cartflip(cartidx::CartesianIndex) = CartesianIndex(cartidx[2], cartidx[1])
        for idx in pairidxs
            if validpairs[idx] != 0
                # Enter twice, once in the upper triangle and once in the lower
                validpairs[idx] = 1
                validpairs[cartflip(idx)] = 1
            end
        end

        # 0 out anything that has not been seen yet and was not among the valid pairs
        netpairs[netpairs .== -1] .= 0
    end

    # Any pairs that still have not been seen are valid
    validpairs[validpairs .== -1] .= 1

    return validpairs
end
findvalidpairs(net::HybridNetwork, namelist::AbstractVector{<:AbstractString}) = findvalidpairs([net], namelist)


"""
    findsiblingpairs(net::HybridNetwork)

Finds sibling pairings for all the leaves in a single network.
These pairs are valid for `net` but may not be valid when the
    other constraint networks are also considered.
Returns a vector of tuples of nodes corresponding to siblings.
"""
function findsiblingpairs(net::HybridNetwork)
    pairs = Vector{Tuple{Node, Node}}()

    # new, simpler method just using Graphs.jl
    hybedges = [edge for edge in net.edge if (edge.hybrid && !edge.isMajor)]
    if length(hybedges) == 0 hybedges = [nothing] end
    for (nodei, nodej) in combinations(net.leaf, 2)
        for hybedge in hybedges
            graph = Graph(net, includeminoredges=false, alwaysinclude=hybedge)
            removeredundantedges!(graph, keeproot=net)

            idxnodei = findfirst(net.node .== [nodei])
            idxnodej = findfirst(net.node .== [nodej])
            edgepath = a_star(graph, idxnodei, idxnodej)

            nodesinpath = Array{Node}(undef, length(edgepath)+1)
            for (i, gedge) in enumerate(edgepath)
                srcnode = net.node[gedge.src]
                dstnode = net.node[gedge.dst]

                if i == 1 nodesinpath[1] = net.node[gedge.src] end
                nodesinpath[i+1] = dstnode
            end

            if length(edgepath) == 2 || length(edgepath) == 1 || (length(edgepath) == 3 && any([(hybedge in n.edge) for n in nodesinpath]))
                push!(pairs, (nodei, nodej))
                break
            end
        end
    end

    return pairs
end


"""

Helper function for manipulating graph adjacency lists.
"""
function vecreplace!(vec::Vector{T}, pattern::T, with::T) where T
    for (i, val) in enumerate(vec)
        if val == pattern
            vec[i] = with
        end
    end
end


"""
A helper function that finds other nodes in the net that `leaf` belongs to that may be its sibling.
This function sees past reticulations, so that e.g. in ((A,(B,#H1)),((C)#H1,D)) B and C will be
    marked as siblings.
"""
function getsiblingcandidates(leaf::Node)
    # NOTE: `children` is used in a loose, unrooted sense here to mean a
    #       node connected to the `parent` by an edge

    # one level up
    parentnode = getnodes(leaf)[1]
    children = getnodes(parentnode)
    followedhyb = repeat([false], length(children))
    visited = Set{Node}([parentnode])
    
    # Pre-process the `children` list so that we can see past reticulations
    keeplooping = true
    while keeplooping
        keeplooping = false
        for (i, child) in enumerate(children)
            if child in visited continue end
            push!(visited, child)
            nextchildren = getnodes(child)

            # If this case is met, the node `child` is an internal node that defined
            # a reticulation. In this case, sibling relationships can still exist 1
            # level deeper, so we iterate one level deeper.
            ishybrid = [nextchild.hybrid || any([e.hybrid for e in nextchild.edge]) for nextchild in nextchildren]
            isunvisited = .![nc in visited for nc in nextchildren]

            if child.hybrid || (any(ishybrid) & any(isunvisited))
                # If one of the children of `child` is a hybrid, then `child` is an internal
                # node corresponding to a reticulation, and we want to be able to see siblings
                # through it. So, we replace it with any of its children that are NOT hybrids.
                # NOT b/c we don't want to follow a potentially large path of reticulations,
                # we just want to maintain sibling status for taxa that would be siblings if
                # the network were pruned of its reticulations.
                keeplooping = true  # go through at least one more time
                
                idx = isunvisited
                if followedhyb[i]
                    idx .&= .![e.hybrid && !e.isMajor for e in child.edge]
                end

                for (j, nextchild) in enumerate(nextchildren[idx])
                    # If we had to cross a reticulation to get to this child,
                    # then we mark it as visited so that we don't continue
                    # our search in that direction anymore
                    if j == 1
                        children[i] = nextchild
                        followedhyb[i] = true
                    else
                        push!(children, nextchild)
                        push!(followedhyb, true)

                        if followedhyb[i]
                            push!(visited, nextchild)
                        end
                    end
                end
            elseif followedhyb[i]
                for n in nextchildren
                    push!(children, n)
                    push!(visited, n)
                end
            end
        end
    end

    return children
end

getnodes(n::Node) = reduce(vcat, [[child for child in e.node if child != n] for e in n.edge], init=Vector{Node}())