# New approach: have tuples with (HybridNetwork, Vector{LineageNode})
#               so that instead of making fundamental changes to
#               HybridNetwork, Vector{LineageNode} is just a 1-to-1
#               vector mapping onto net.node


function getparentaltrees(newick::AbstractString; safe::Bool=true)
    return getparentaltrees(readTopology(newick))
end


# Gets the parental trees for the given network.
# !! TREES MAY NOT BE UNIQUE !!
function getparentaltrees(net::HybridNetwork; safe::Bool=true)
    if safe net = deepcopy(net) end
    _overwritemissinggammas!(net)
    _overwritemissingbranchlengths!(net)

    # Save original leaf labels
    leaflabs = [leaf.name for leaf in net.leaf]

    # NETWORK PRE-PROCESSING
    ldict = LDict()   # nodes & nets will be deep copied, 
                            # so this never needs to be deep copied
    i = 1
    for node in net.node
        if node.leaf
            ldict[node] = LineageNode(i)
            node.name = string(i)
            i += 1
        else
            ldict[node] = nothing
        end
    end

    ipt = IPT(net)
    ptrees = _getparentaltrees(ipt, ldict)

    # Expand the parental trees from `ldict` so that the `HybridNetwork`
    # topology reflects the actual parental tree topology
    ptrees = [_expandIPT(pt, ldict) for pt in ptrees]

    # TODO: Condense non-unique parental trees

    # Put the original leaf labels back on the trees
    for pt in ptrees
        for leaf in top(pt).leaf
            leaf.name = leaflabs[parse(Int64, leaf.name)]
        end
    end

    return ptrees, ldict
end

function _getparentaltrees(initnet::IPT, lineagedict::LDict)
    workingset = Queue{IPT}()
    enqueue!(workingset, initnet)
    parentaltrees = Vector{IPT}()

    while !isempty(workingset)
        ipt = dequeue!(workingset)
        net = top(ipt)

        # 0. No more hybrids? This net is done!
        if length(net.hybrid) == 0
            push!(parentaltrees, ipt)
            continue
        end

        # 1. Find a hybrid node with no other hybrids below it
        node, nodeidx = _getnexthybrid(net)
        divisions = Vector{IPT}()

        # 2a. If there is already a LineageNode associated with this node,
        #     then we are ready to condition on the reticulation
        if lineagedict[node] !== nothing
            print("Conditioning on reticulations - ")
            divisions = _conditiononreticulation(ipt, node, lineagedict)
            println("done.")

        # 2b. If there is not a LineageNode associated with thisnode,
        #     we must first condition on the coalescent events
        else
            print("Conditioning on coalescences - ")
            divisions = _conditiononcoalescences(ipt, node, lineagedict)
            println("done.")
        end

        # 4. Repeat with the resultant networks
        for d in divisions enqueue!(workingset, d) end
    end

    return parentaltrees
end

# Overwrites any missing/unspecified gammas for hybrids in `net`
function _overwritemissinggammas!(net::HybridNetwork)
    invalidgamma(edge::Edge) = edge.gamma == -1. || edge.gamma > 1. || edge.gamma < 0.
    
    warn = false
    for hyb in net.hybrid
        major = getparentedge(hyb)
        minor = getparentedgeminor(hyb)
        if invalidgamma(major)
            warn = true
            if invalidgamma(minor)
                major.gamma = 0.5
                minor.gamma = 0.5
            else
                major.gamma = 1 - minor.gamma
            end
        elseif invalidgamma(minor)
            warn = true
            minor.gamma = 1 - major.gamma
        end
    end

    if warn @warn "One or more invalid gammas found and replaced in original topology." end
end

function _overwritemissingbranchlengths!(net::HybridNetwork; replacement::Real=1.)
    warn = false
    for edge in net.edge
        if edge.length < 0.
            warn = true
            edge.length = replacement
        end
    end

    if warn @warn "One or more branch lengths were unspecified and have been replaced with 1.0." end
end

# Takes an incomplete parental tree and expands its stored LineageNode attributes
# into a complete parental tree (returned as a `HybridNetwork` object)
# TODO: implement this
function _expandIPT(tree::IPT, ldict::LDict)
    return tree
end

# DOES NOT WORK AS WRITTEN
# Problem: investigates `divisions` ONLY at each node in `nodes`
#          if 2 divisions are identical **at this node** but distinct elsewhere,
#          they will be merged here
function _condensedivisions(divisions::Vector{HybridNetwork}, nodes::Vector{Node}, ldict::LDict)
    println("\nCondensing")
    if length(divisions) < 2 return divisions end
    println("Past length 1 check")

    # Compare the LineageNode entries for each division at `node`
    # Merges divisions w/ unique `ldict[node]` entries
    i = 1
    j = 2
    while true
        # print("(", i, ",", j, "): ")
        # print(ldict[nodes[i]])
        # print(" == ")
        # println(ldict[nodes[j]])
        if ldict[nodes[i]] == ldict[nodes[j]]
            # MERGE PROBABILITIES HERE
            deleteat!(nodes, j)
            deleteat!(divisions, j)
            # println("Deleted")
        else
            j += 1
        end

        if j > length(divisions)
            i += 1
            j = i + 1
            if i >= length(divisions) break end
        end
    end

    return divisions
end

# Conditions the given network on the given reticulation
function _conditiononreticulation(ipt::IPT, hyb::Node, lineagedict::LDict)
    # If everything is moving as expected, the children should be exclusively
    # one leaf node (signifying a single lineage) or one LineageNode
    net = top(ipt)
    child = lineagedict[hyb]   # this will throw an error if more than 1 child exists
    hybidx = findfirst(net.node .== [hyb])

    retlist = Vector{IPT}()
    for set in powerset([1:nlineages(child);])
        # `set` contains the indices of the children that go up the "left" reticulation
        # while its `setdiff` is those the go up the "right" reticulation
        majorline = LineageNode(lineages(child)[set])
        minorline = LineageNode(lineages(child)[setdiff(1:nlineages(child), set)])

        # Probability of the given split:
        #   γ: left inheritance probability
        #   nmaj: number of lineages going left
        #   nmin: number of lineages going right
        #   P(split) = γ^nmaj * (1-γ)^nmin
        γ = getparentedge(hyb).gamma
        nmin = length(set)
        tot = nlineages(child)  # nmaj = tot-nmin
        splitprob = γ^(tot-nmin) * (1-γ)^nmin
        newprob = prob(ipt) * splitprob

        newnet = deepcopy(net)
        newhyb = newnet.node[hybidx]
        copyldictcontents!(net, newnet, lineagedict)

        # Split the reticulation into 2 tree-like splits
        splitreticulation!(newnet, newhyb, hybidx, majorline, minorline, lineagedict)
        push!(retlist, IPT(newnet, newprob))
    end

    return retlist
end

# Conditions the given network on the various coalescent possibilities below `node`
function _conditiononcoalescences(ipt::IPT, node::Node, lineagedict::LDict)
    net = top(ipt)
    if lineagedict[node] !== nothing return ipt end    # this means we'll never hit leaves
    
    children = getchildren(node)
    length(children) <= 2 || error("Polytomies are not accounted for. Found "*string(length(children))*" children of node #"*string(node.number))

    # Check that the children are good; if not, recurse
    ipts = Vector{IPT}([ipt])
    if lineagedict[children[1]] === nothing
        ipts = reduce(vcat, [_conditiononcoalescences(ipt, children[1], lineagedict) for ipt in ipts])
    end
    if length(children) > 1 && lineagedict[children[2]] === nothing
        ipts = reduce(vcat, [_conditiononcoalescences(ipt, children[2], lineagedict) for ipt in ipts])
    end

    # In the simplest case, this only has 1 net and `node` is above 2 leaves
    retlist = Vector{IPT}()
    # markednodes = Vector{Node}()    # used for condensing later
    for tempipt in ipts
        tempnet = top(tempipt)

        # when we copy nets, we need to somehow make sure that lower down references are not lost
        tempnodeidx = findfirst([n.number == node.number for n in tempnet.node])
        tempnode = tempnet.node[tempnodeidx]
        children = getchildren(tempnode)

        # We are thinking about each of the children coming up from their node and into this node.
        # From that perspective, all the lineages in each of the children nodes have _separate_
        # chances to merge while moving up their branch. So, the resultant node `tempnode`
        # potentially has any pairing of `opts1` and `opts2`
        opts1, probs1 = getcoalescentcombos(lineagedict[children[1]], getparentedge(children[1]).length)
        if length(children) > 1
            # println("children[2] address: "*repr(UInt64(pointer_from_objref(children[2]))))
            opts2, probs2 = getcoalescentcombos(lineagedict[children[2]], getparentedge(children[2]).length)
        else
            opts2 = [LineageNode()]
            probs2 = [1]
        end
        
        for (i, (opti, probi)) in enumerate(zip(opts1, probs1))
            for (j, (optj, probj)) in enumerate(zip(opts2, probs2))
                newnet = deepcopy(tempnet)
                copyldictcontents!(tempnet, newnet, lineagedict)
                newprob = prob(tempipt) * probi * probj

                lineagedict[newnet.node[tempnodeidx]] = LineageNode(opti, optj)
                push!(retlist, IPT(newnet, newprob))
                # push!(markednodes, newnet.node[tempnodeidx])
            end
        end
    end

    # # Condense down to unique parental trees
    # _condensedivisions(retlist, markednodes, lineagedict)
    return retlist
end

# Gets all the _deep_ children of `hyb`. I.e., not the nodes immediately proceeding it,
# but all of the furthest tips that it can reach (not including any intermediate nodes)
function _getchildrendeep(hyb::Node)
    children = getchildren(hyb)
    returnlist = Vector{Node}([child for child in children if typeof(child) <: LineageNode || child.leaf])
    for child in setdiff(children, returnlist)
        returnlist = vcat(returnlist, _getchildrendeep(child))
    end

    return returnlist
end

# Finds a hybrid node with no other hybrid nodes below it
function _getnexthybrid(net::HybridNetwork)
    for (i, hyb) in enumerate(net.hybrid)
        if _noreticsbelow(hyb)
            return hyb, i
        end
    end

    return nothing
end

# Checks whether there are any hybrid nodes below `hyb` in the network
function _noreticsbelow(node::Node)
    if node.leaf return true end

    children = getchildren(node)
    
    # Search breadth first, not depth first
    for child in children
        if child.hybrid return false end
    end

    for child in children
        if !_noreticsbelow(child) return false end
    end

    return true
end