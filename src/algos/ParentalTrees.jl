# Good example net to work with: ((A,((B,(C,D)))#H1),(#H1,E));

using PhyloNetworks
import Combinators: powerset
import DataStructures: Queue


function getParentalTrees(net::HybridNetwork)
    # NETWORK PRE-PROCESSING
    # 

    return _getParentalTrees(deepcopy(net))
end

function _getParentalTrees(net::HybridNetwork)
    while length(net.hybrid) != 0
        # 1. Find a hybrid node with no other hybrids below it
        node = _getNextHybrid(net)
        divisions = Vector{HybridNetwork}()
    
        # 2a. If there is only 1 LineageNode below the reticulation,
        #    condition on the reticulation
        children = getchildren(node)
        if length(children) == 1 && (typeof(children) <: LineageNode || children[1].leaf)
            divisions = _conditionOnReticulation(net, node)

        # 2b. If there are multiple LineageNodes below the reticulation,
        #    condition on coalescent events
        else
            divisions = _conditionOnCoalescences(net, node)
        end

        # 3. Repeat with the resultant networks
    end

    return net
end

# Conditions the given network on the given reticulation
function _conditionOnReticulation(net::HybridNetwork, hyb::PhyloNetworks.ANode)
    # If everything is moving as expected, the children should be exclusively
    # one leaf node (signifying a single lineage) or one LineageNode
    child = getchild(hyb)   # this will throw an error if more than 1 child exists

    if typeof(child) <: PhyloNetworks.Node || nlineages(child) == 1
        # Only 2 outcomes: follow the left or follow the right

    else
        typeof(child) <: LineageNode || error("Expected a LineageNode, instead found a "*str(typeof(child)))
        # 2^[nlineages(child)] outcomes from here
        # Need to split among the set of all left/right combinations for each lineage

        for set in powerset([1:nlineages(child);])
            # `set` contains the indices of the children that go up the "left" reticulation
            # while its `setdiff` is those the go up the "right" reticulation
            leftLineages = lineages(child)[set]
            rightLineages = lineages(child)[setdiff(1:nlineages(child), set)]
        end
    end
end

function _conditionOnCoalescences(net::HybridNetwork, node::PhyloNetworks.ANode)
    return _conditionOnCoalescences([net], node)
end

# Conditions the given network on the various coalescent possibilities below `node`
function _conditionOnCoalescences(nets::Vector{HybridNetwork}, node::PhyloNetworks.ANode)
    children = getchildren(node)
    length(children) == 2 || error("Polytomies are not accounted for.")

    if !(typeof(children[1]) <: LineageNode || children[1].leaf)
        nets = reduce(vcat, [_conditionOnCoalescences(net, children[1]) for net in nets])  # ISSUE: `node` will (I think) be different in each of these networks,
                                                                                           # so may have to store the data as a vector of tuples (net, node)
    end

    if !(typeof(children[2]) <: LineageNode || children[2].leaf)
        nets = reduce(vcat, [_conditionOnCoalescences(net, children[2]) for net in nets])
    end

    # Now the left and right nodes should both either be leaves or LineageNodes
    if typeof(children[1]) <: LineageNode
        if typeof(chidren[2]) <: LineageNode
            # both are LineageNodes

        else
            # children[1] is LineageNode, children[2] is leaf

        end
    else
        if typeof(children[2]) <: LineageNode
            # children[1] is leaf, children[2] is LineageNode

        else
            # easiest case; both are leaf nodes

        end
    end
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