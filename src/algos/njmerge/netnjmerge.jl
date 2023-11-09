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
    end

    # Empty network
    subnets = [SubNet(i, names) for i in 1:n]
    edgenum = 0
    reticnum = 0
    reticmap = Vector{Tuple{Edge, Edge}}()

    # Keeping track of the algorithm
    possible_siblings = findvalidpairs(constraints, names)

    # Used to efficiently compute Q
    Dsums = sum(D, dims=1)

    while n > 2
        # Find optimal (i, j) idx pair for matrix Q
        i, j = findoptQidx(D, idxpairs)

        # connect subnets i and j
        # TODO: before implementing this section, sketch out
        #       what it should look like to connect 2 subnets
        #       (review nj! code)
        subnets[i] = mergesubnets!(subnets[i], subnets[j])

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

    finalnet = mergesubnets!(subnets[i], subnets[j])
    finalnet = HybridNetwork(finalnet)

    # Place the reticulations we've been keeping track of
    placeretics!(finalnet, reticmap)

    return finalnet
end


"""
    findoptQ(D::Matrix{Float64}, idxpairs::Vector{Tuple{<:Integer, <:Integer}})

Finds the minimizer (i*, j*) among all pairs (i, j) in idxpairs for Q, a matrix computed from D.
"""
function findoptQidx(D::Matrix{Float64}, idxpairs::Vector{Tuple{<:Integer, <:Integer}})
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
            ishybrid = [nextchild.hybrid for nextchild in nextchildren]
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
                    idx .&= .![e.hybrid for e in child.edge]
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