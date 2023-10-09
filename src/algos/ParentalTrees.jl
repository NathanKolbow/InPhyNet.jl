function getparentaltrees(newick::AbstractString; safe::Bool=true)
    return getparentaltrees(readTopology(newick))
end

"""
    getparentaltrees(net::HybridNetwork; safe::Bool=true)

Gets the set of parental trees for network `net` and their associated
probabilities under the multispecies network coalescent (MSNC) model.

Trees may not be unique.
"""
function getparentaltrees(net::HybridNetwork; safe::Bool=true)
    if safe net = deepcopy(net) end

    ipt, ldict, labelmap = _initparentaltreealgo(net)
    ptrees, complog = _getparentaltrees(ipt, ldict, usecomplog=true)

    # Expand the parental trees from `ldict` so that the `HybridNetwork`
    # topology reflects the actual parental tree topology
    ptrees = [_expandIPT(pt, ldict) for pt in ptrees]

    # TODO: Condense non-unique parental trees

    # Put the original leaf labels back on the trees
    for pt in ptrees
        for leaf in top(pt).leaf
            leaf.name = labelmap[parse(Int64, leaf.name)]
        end
    end

    return ptrees, ldict
end

"""
    _initparentaltreealgo(net::HybridNetwork)

Takes the given network and returns initialization values for the parental tree algorithm.

# Return Values
1. InterimParentalTree generated from the given network
2. Initialized LDict
"""
function _initparentaltreealgo(net::HybridNetwork; neverwarn=true)
    # Fill in missing gammas and branch lengths
    _overwritemissinggammas!(net, neverwarn=neverwarn)
    _overwritemissingbranchlengths!(net, neverwarn=neverwarn)

    # Fuse redundant edges
    _fuseredundantedges!(net)

    # Save original leaf labels
    labelmap = Dict{Int64, String}()

    # NETWORK PRE-PROCESSING
    ldict = LDict()   # nodes & nets will be deep copied, 
                      # so this never needs to be deep copied
    i = 1
    for node in net.node
        if node.leaf
            # Save the tip label
            labelmap[i] = node.name

            # Initialize the LDict
            ldict[node] = LineageNode(i)

            # Overwrite the name (not really necessary)
            node.name = string(i)
            # Increment
            i += 1
        else
            # Initialize the LDict entry
            ldict[node] = nothing
        end
    end

    return IPT(net), ldict, labelmap
end

"""
    _fuseredundantedges!(net::HybridNetwork)

Traverses `net` for redundant nodes & edges. E.g. (A:1,(B:1):1) is
redundant and (A:1,B:2) makes things computationally more efficient.
Used when initializing the parental tree algorithm. Care is taken
not to create more redundant nodes & edges throughout the algorithm.
"""
function _fuseredundantedges!(net::HybridNetwork)
    return nothing
end

function _getparentaltrees(initnet::IPT, ldict::LDict; usecomplog=true)
    workingset = Queue{IPT}()
    enqueue!(workingset, initnet)
    parentaltrees = Vector{IPT}()
    complog = ifelse(usecomplog, CompDict(), nothing)   # for re-using computations

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
        if ldict[node] !== nothing
            divisions = _conditiononreticulation(ipt, node, ldict)

            newsum = sum([prob(t) for t in divisions])
            oldprob = prob(ipt)
            abs(newsum - oldprob) < 1e-16 || error("Probs don't line up after splitting retics; newsum="*string(newsum)*", oldprob="*string(oldprob))
        # 2b. If there is not a LineageNode associated with thisnode,
        #     we must first condition on the coalescent events
        else
            divisions = _conditiononcoalescences(ipt, node, ldict, complog=complog)

            newsum = sum([prob(t) for t in divisions])
            oldprob = prob(ipt)
            abs(newsum - oldprob) < 1e-16 || error("Probs don't line up after coalescing from node"*string(node.name)*"; newsum="*string(newsum)*", oldprob="*string(oldprob)*". length(divisions): "*string(length(divisions)))
        end

        # 4. Repeat with the resultant networks
        for d in divisions enqueue!(workingset, d) end
    end

    return parentaltrees, complog
end

# Overwrites any missing/unspecified gammas for hybrids in `net`
function _overwritemissinggammas!(net::HybridNetwork; neverwarn::Bool=true)
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

    if !neverwarn && warn @warn "One or more invalid gammas found and replaced in original topology." end
end

function _overwritemissingbranchlengths!(net::HybridNetwork; replacement::Real=1., neverwarn::Bool=true)
    warn = false
    for edge in net.edge
        if edge.length < 0.
            warn = true
            edge.length = replacement
        end
    end

    if !neverwarn && warn @warn "One or more branch lengths were unspecified and have been replaced with 1.0." end
end

# Takes an incomplete parental tree and expands its stored LineageNode attributes
# into a complete parental tree (returned as a `HybridNetwork` object)
# TODO: implement this
function _expandIPT(tree::IPT, ldict::LDict)
    return tree
end

# Used to use this, but was implemented improperly.
# Awaiting reimplementation.
# Arguments need not remain the same.
function _condensedivisions(something)
    return nothing
end

# Conditions the given network on the given reticulation
function _conditiononreticulation(ipt::IPT, hyb::Node, ldict::LDict)
    # If everything is moving as expected, the children should be exclusively
    # one leaf node (signifying a single lineage) or one LineageNode
    net = top(ipt)
    child = ldict[hyb]   # this will throw an error if more than 1 child exists
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
        copyldictcontents!(net, newnet, ldict)

        # Split the reticulation into 2 tree-like splits
        splitreticulation!(newnet, newhyb, hybidx, majorline, minorline, ldict)
        push!(retlist, IPT(newnet, newprob))
    end

    return retlist
end

# Conditions the given network on the various coalescent possibilities below `node`
function _conditiononcoalescences(ipt::IPT, node::Node, ldict::LDict; complog::Union{CompDict,Nothing}=nothing)
    net = top(ipt)
    if ldict[node] !== nothing return ipt end    # this means we'll never hit leaves
    
    children = getchildren(node)
    length(children) <= 2 || error("Polytomies are not accounted for. Found "*string(length(children))*" children of node #"*string(node.number))

    # Check that the children are good; if not, recurse
    ipts = Vector{IPT}([ipt])
    if ldict[children[1]] === nothing
        ipts = reduce(vcat, [_conditiononcoalescences(ipt, children[1], ldict, complog=complog) for ipt in ipts])
    end
    if length(children) > 1 && ldict[children[2]] === nothing
        ipts = reduce(vcat, [_conditiononcoalescences(ipt, children[2], ldict, complog=complog) for ipt in ipts])
    end

    # In the simplest case, this only has 1 net and `node` is above 2 leaves
    retlist = Vector{IPT}()
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
        opts1, probs1 = getcoalescentcombos(ldict[children[1]], getparentedge(children[1]).length, complog=complog)
        if length(children) > 1
            opts2, probs2 = getcoalescentcombos(ldict[children[2]], getparentedge(children[2]).length, complog=complog)
        else
            opts2 = [LineageNode()]
            probs2 = [BigFloat(1.)]
        end

        for (i, (opti, probi)) in enumerate(zip(opts1, probs1))
            for (j, (optj, probj)) in enumerate(zip(opts2, probs2))
                newnet = deepcopy(tempnet)
                copyldictcontents!(tempnet, newnet, ldict)
                newprob = prob(tempipt) * probi * probj

                ldict[newnet.node[tempnodeidx]] = LineageNode(opti, optj)
                for child in getchildren(newnet.node[tempnodeidx]) delete!(ldict, child) end
                push!(retlist, IPT(newnet, newprob))
            end
        end
    end

    # TODO: condense the retlist right here! everything in `retlist` came from the same initial `ipt` so will
    #       be identical up to their differences at `node`. So, we can condense anything that is the same there
    #       AND we can condense the ambiguous coalescences out (e.g. [[1, 2, 3]] w/ prob 0.3 split out to [[1, 2], 3],
    #       [[1, 3], 2], [1, [2, 3]] giving probability 0.1 to each)

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
    if length(net.hybrid) == 1 return net.hybrid[1], 1 end
    for (i, hyb) in enumerate(net.hybrid)
        if _noreticsbelow(hyb)
            return hyb, i
        end
    end

    return nothing
end
_getnexthybrid(ipt::IPT) = _getnexthybrid(top(ipt))

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