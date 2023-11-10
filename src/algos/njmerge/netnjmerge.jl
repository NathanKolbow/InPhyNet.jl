import PhyloNetworks: deleteNode!, deleteEdge!

# Source for the network version of NJ merge
function netnj!(D::Matrix{Float64}, constraints::Vector{HybridNetwork};
    names::AbstractVector{<:AbstractString}=String[])
    
    PhyloNetworks.check_distance_matrix(D)
    n = size(D, 1)

    # If names are not provided
    if isempty(names)
        names = string.(1:n)
    end
    if length(names) != n
        m = length(names)
        throw(ArgumentError("D has dimensions $n x $n but only $m names were provided."))
    elseif length(unique(names)) != length(names)
        throw(ArgumentError("Names must be unique."))
    end

    # Empty network
    subnets = Vector{SubNet}([SubNet(i, names[i]) for i in 1:n])
    reticmap = ReticMap(constraints)

    # Used to efficiently compute Q
    Dsums = sum(D, dims=1)

    while n > 2
        possible_siblings = findvalidpairs(constraints, names)
        
        # Find optimal (i, j) idx pair for matrix Q
        i, j = findoptQidx(D, possible_siblings)

        # connect subnets i and j
        # TODO: before implementing this section, sketch out
        #       what it should look like to connect 2 subnets
        #       (review nj! code)
        subnets[i], edgei, edgej = mergesubnets!(subnets[i], subnets[j])
        updateconstraints!(names[i], names[j], constraints, reticmap, edgei, edgej)

        # collapse taxa i into j
        for l in 1:n
            if l != i && l != j
                D[l, i] = D[i, l] = (D[l, i] + D[j, l] - D[i, j]) / 2
            end
        end

        # Remove data elements that corresponded to `j`
        idxfilter = [1:(j-1); (j+1):n]
        D = view(D, idxfilter, idxfilter)   # remove j from D
        subnets = view(subnets, idxfilter)

        n -= 1
    end

    return mergesubnets!(subnets[1], subnets[2])[1]

    finalnet, _, _ = mergesubnets!(subnets[1], subnets[2])
    finalnet = HybridNetwork(finalnet)

    # Place the reticulations we've been keeping track of
    # placeretics!(finalnet, reticmap)

    return finalnet
end


"""

Updates `net` to include the reticulations that we've kept track of along the way
in our algo but haven't placed yet.
"""
function placeretics!(net::HybridNetwork, reticmap::ReticMap)
    error("not yet implemented")
end


"""

Updates constraint networks after (i, j) with names (nodenamei, nodenamej) have been
merged. Also updates `reticmap` to keep track of any reticulations that get removed
in this process.

By convention we keep `nodenamei` and replace node names with `nodenamej`
"""
function updateconstraints!(nodenamei::AbstractString, nodenamej::AbstractString, 
    constraints::Vector{HybridNetwork}, reticmap::ReticMap, subnetedgei::Edge, subnetedgej::Edge)

    for net in constraints
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
    length(parentsi) == 1 || error("Found >1 nodes above a leaf?")    # sanity check, remove when finalized
    length(parentsj) == 1 || error("Found >1 nodes above a leaf?")    # sanity check; remove when finalized
    parentsi = parentsi[1]
    parentsj = parentsj[1]

    if parentsi == parentsj
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
        error("not tested yet")
        error("this is gonna need a good amount of testing...")

        # reticulations: make our way up and log retics in `reticmap`

        # just traverse away from nodes i and j, marking nodes as visited as we go,
        # and when we come to a node w/ 2 parents taking the `isMajor` edge, otherwise
        # the node will have 1 parent and we can just continue onwards
        prevnodei = nothing
        curri = nodei
        prevnodej = nothing
        currj = nodej

        visitedbyi = Vector{Node}([nodei])  # vectors instead of sets b/c order matters
        visitedbyj = Vector{Node}([nodej])  # vectors instead of sets b/c order matters
        edgepathi = Vector{Edge}()
        edgepathj = Vector{Edge}()
        newtip = nodei

        while true
            # Move to next node
            curri, fromedgei = getnextmajornode(nodei, prevnodei, reticmap, subnetedgei)
            currj, fromedgej = getnextmajornode(nodej, prevnodej, reticmap, subnetedgej)
            push!(edgepathi, fromedgei)
            push!(edgepathj, fromedgej)

            # Add new nodes to visited sets
            push!(visitedbyi, nodei)
            push!(visitedbyj, nodej)

            # If we've crossed paths, break
            if curri in visitedbyj
                newtip = curri
                break
            elseif currj in visitedbyi
                newtip = currj
                break
            end
        end

        # now sever all the nodes we've passed through to this point from the network
        tipidxini = findfirst(visitedbyi .== [newtip])
        tipidxinj = findfirst(visitedbyj .== [newtip])
        purgenodes = union(visitedbyi[1:(tipidxini-1)], visitedbyj[1:(tipidxinj-1)])
        purgeedges = union(edgepathi[1:tipidxini], edgepathj[1:tipidxinj])

        for node in purgenodes
            deleteNode!(net, node)
        end
        for edge in purgeedges
            for node in edge.node
                node.edge = filter(e -> e != edge, node.edge)
            end
            deleteEdge!(net, edge)
        end

        newtip.leaf = true
        newtip.name = nodei.name
    end
end


function getnextmajornode(currnode::Node, prevnode::Node, reticmap::ReticMap, subnetedge::Edge)
    candidateedges = filter(e -> !(prevnode in e.node), node.edge)  # everything that doesn't go backwards
    nextmajoredge = filter(e -> !e.hybrid || e.isMajor, candidateedges)
    length(nextmajoredge) == 1 || error("Should only be 1 edge matching these conditions")  # sanity check

    if nextmajoredge.hybrid
        # then we need to push its minor edge partner to `reticmap`
        nextminoredge = filter(e -> e.hybrid && !e.isMajor, candidateedges)
        logretic(reticmap, nextminoredge, subnetedge)
    end

    return filter(n -> n != currnode, nextmajoredge.node)[1], nextmajoredge
end


"""
    findoptQ(D::Matrix{Float64}, idxpairs::Vector{Tuple{<:Integer, <:Integer}})

Finds the minimizer (i*, j*) among all pairs (i, j) in idxpairs for Q, a matrix computed from D.
"""
function findoptQidx(D::AbstractMatrix{Float64}, validpairs::Matrix{<:Integer})
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
        end
    end
    
    return minidx
end


"""
    findvalidpairs(constraints::Vector{HybridNetwork}, names::AbstractVector{<:AbstractString})

Finds all valid sibling pairs among the constraint networks.
"""
function findvalidpairs(constraints::Vector{HybridNetwork}, names::AbstractVector{<:AbstractString})
    n = length(names)

    # Shorthand functions for name lookups (we'll be doing a lot of these if there are many constraints)
    namedict = Dict{AbstractString, Int64}([name => i for (i, name) in enumerate(names)])
    idx(name::AbstractString) = namedict[name]

    # initialize matrix
    validpairs = Matrix{Int64}(undef, n, n)     # -1: pair not seen together yet
                                                #  0: pair invalid
                                                #  1: pair valid
    validpairs .= -1

    # go through the constraint networks and validate/invalidate pairs
    for net in constraints
        leafidxs = [idx(leaf.name) for leaf in net.leaf]

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
        netpairs = view(validpairs, leafidxs, leafidxs)
        netpairs[netpairs .== -1] .= 0
    end

    # Any pairs that still have not been seen are valid
    validpairs[validpairs .== -1] .= 1

    return validpairs
end
findvalidpairs(net::HybridNetwork, names::AbstractVector{<:AbstractString}) = findvalidpairs([net], names)


"""
    findsiblingpairs(net::HybridNetwork)

Finds sibling pairings for all the leaves in a single network.
These pairs are valid for `net` but may not be valid when the
    other constraint networks are also considered.
Returns a vector of tuples of nodes corresponding to siblings.
"""
function findsiblingpairs(net::HybridNetwork)
    pairs = Vector{Tuple{Node, Node}}()
    already_examined = Array{Node}(undef, net.numTaxa)
    l = 0
    
    for leaf in net.leaf
        # Get the list of potential siblings
        children = getsiblingcandidates(leaf)

        # Now hybrids have already been processed out, so we can treat
        # any `child` in `children` as an actual potential sibling to `leaf`
        for child in children
            # If this child is (1.) also a leaf, (2.) not the node we're currently looking at
            # and (3.) hasn't already been looked at (to avoid duplicates), then add it
            if child.leaf && child != leaf && !(child in already_examined[1:l])
                push!(pairs, (leaf, child))
            end
        end

        # update the list of already examined leaves
        already_examined[l+1] = leaf
        l += 1
    end

    return pairs
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

getnodes(n::Node) = reduce(vcat, [[child for child in e.node if child != n] for e in n.edge])

# Test for `findsiblingpairs`
# nn = readTopology("((A,(B,#H1)),((C)#H1,D));")
# pairs = findsiblingpairs(nn)
# length(pairs) == 3 || error("test failed in netnjmerge.jl")