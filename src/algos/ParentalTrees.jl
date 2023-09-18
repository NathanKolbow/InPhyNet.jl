# New approach: have tuples with (HybridNetwork, Vector{LineageNode})
#               so that instead of making fundamental changes to
#               HybridNetwork, Vector{LineageNode} is just a 1-to-1
#               vector mapping onto net.node

using PhyloNetworks
import Combinatorics: powerset
using DataStructures


function getParentalTrees(net::HybridNetwork; safe::Bool=true)
    if safe net = deepcopy(net) end

    # NETWORK PRE-PROCESSING
    ldict = LDict()   # nodes & nets will be deep copied, 
                            # so this never needs to be deep copied
    i = 1
    for node in net.node
        if node.leaf
            ldict[node] = LineageNode(i)
            i += 1
        else
            ldict[node] = nothing
        end
    end

    return _getParentalTrees(net, ldict)
end

function _getParentalTrees(initnet::HybridNetwork, lineagedict::LDict)
    workingset = Queue{HybridNetwork}()
    enqueue!(workingset, initnet)
    parentaltrees = Vector{HybridNetwork}()

    while !isempty(workingset)
        net = dequeue!(workingset)

        # 0. No more hybrids? This net is done!
        if length(net.hybrid) == 0
            push!(parentaltrees, net)
            continue
        end

        # 1. Find a hybrid node with no other hybrids below it
        node = _getNextHybrid(net)
        divisions = Vector{HybridNetwork}()

        # 2a. If there is already a LineageNode associated with this node,
        #     then we are ready to condition on the reticulation
        if lineagedict[node] !== nothing
            print("Conditioning on reticulations - ")
            divisions = _conditionOnReticulation(net, node, lineagedict)
            println("done.")

        # 2b. If there is not a LineageNode associated with thisnode,
        #     we must first condition on the coalescent events
        else
            print("Conditioning on coalescences - ")
            divisions = _conditionOnCoalescences(net, node, lineagedict)
            println("done.")
        end

        # 3. Repeat with the resultant networks
        for d in divisions enqueue!(workingset, d) end
    end

    return parentaltrees
end

# Conditions the given network on the given reticulation
function _conditionOnReticulation(net::HybridNetwork, hyb::PhyloNetworks.ANode, lineagedict::LDict)
    # If everything is moving as expected, the children should be exclusively
    # one leaf node (signifying a single lineage) or one LineageNode
    child = lineagedict[getchild(hyb)]   # this will throw an error if more than 1 child exists
    hybidx = findfirst(net.node .== [hyb])

    retlist = Vector{HybridNetwork}()
    for set in powerset([1:nlineages(child);])
        # `set` contains the indices of the children that go up the "left" reticulation
        # while its `setdiff` is those the go up the "right" reticulation
        leftline = LineageNode(lineages(child)[set])
        rightline = LineageNode(lineages(child)[setdiff(1:nlineages(child), set)])

        newnet = deepcopy(net)
        newhyb = newnet.node[hybidx]
        copyldictcontents!(net, newnet, lineagedict)

        # Split the reticulation into 2 tree-like splits
        splitreticulation!(newnet, newhyb, leftline, rightline, lineagedict)
        push!(retlist, newnet)
    end
    return retlist
end

# Conditions the given network on the various coalescent possibilities below `node`
function _conditionOnCoalescences(net::HybridNetwork, node::PhyloNetworks.ANode, lineagedict::LDict)
    if lineagedict[node] !== nothing return net end    # this means we'll never hit leaves
    
    children = getchildren(node)
    length(children) <= 2 || error("Polytomies are not accounted for. Found "*string(length(children))*" children of node #"*string(node.number))

    # Check that the children are good; if not, recurse
    nets = Vector{HybridNetwork}([net])
    if lineagedict[children[1]] === nothing
        nets = reduce(vcat, [_conditionOnCoalescences(net, children[1], lineagedict) for net in nets])
    end
    if length(children) > 1 && lineagedict[children[2]] === nothing
        nets = reduce(vcat, [_conditionOnCoalescences(net, children[2], lineagedict) for net in nets])
    end

    # In the simplest case, this only has 1 net and `node` is above 2 leaves
    retlist = Vector{HybridNetwork}()
    for tempnet in nets
        # when we copy nets, we need to somehow make sure that lower down references are not lost
        tempnodeidx = findfirst([n.number == node.number for n in tempnet.node])
        tempnode = tempnet.node[tempnodeidx]
        children = getchildren(tempnode)

        # We are thinking about each of the children coming up from their node and into this node.
        # From that perspective, all the lineages in each of the children nodes have _separate_
        # chances to merge while moving up their branch. So, the resultant node `tempnode`
        # potentially has any pairing of `opts1` and `opts2`
        opts1 = getcoalescentcombos(lineagedict[children[1]])
        if length(children) > 1
            opts2 = getcoalescentcombos(lineagedict[children[2]])
        else
            opts2 = [LineageNode()]
        end

        
        for (i, opti) in enumerate(opts1)
            for (j, optj) in enumerate(opts2)
                newnet = deepcopy(tempnet)
                copyldictcontents!(tempnet, newnet, lineagedict)

                lineagedict[newnet.node[tempnodeidx]] = LineageNode(opti, optj)
                push!(retlist, newnet)
            end
        end
    end
    return retlist
end

# Gets all the _deep_ children of `hyb`. I.e., not the nodes immediately proceeding it,
# but all of the furthest tips that it can reach (not including any intermediate nodes)
function _getChildrenDeep(hyb::PhyloNetworks.ANode)
    children = getchildren(hyb)
    returnlist = Vector{PhyloNetworks.ANode}([child for child in children if typeof(child) <: LineageNode || child.leaf])
    for child in setdiff(children, returnlist)
        returnlist = vcat(returnlist, _getChildrenDeep(child))
    end

    return returnlist
end

# Finds a hybrid node with no other hybrid nodes below it
function _getNextHybrid(net::HybridNetwork)
    for hyb in net.hybrid
        if _noReticsBelow(hyb)
            return hyb
        end
    end

    return nothing
end

# Checks whether there are any hybrid nodes below `hyb` in the network
function _noReticsBelow(node::PhyloNetworks.ANode)
    if node.leaf return true end

    children = getchildren(node)
    
    # Search breadth first, not depth first
    for child in children
        if child.hybrid return false end
    end

    for child in children
        if !_noReticsBelow(child) return false end
    end

    return true
end