using PhyloNetworks
import Combinators: powerset



function getParentalTrees(net::HybridNetwork)
    # NETWORK PRE-PROCESSING
    # 

    return _getParentalTrees(net)
end

function _getParentalTrees(net::HybridNetwork)
    while length(net.hybrid) != 0
        # 1. Find a hybrid node with no other hybrids below it
        node = _getNextHybrid(net)
    
        # 2a. If there is only 1 LineageNode below the reticulation,
        #    condition on the reticulation
        if length(getchildren(node)) == 1
            divisions = _conditionOnReticulation(net, node)
        end

        # 2b. If there are multiple LineageNodes below the reticulation,
        #    condition on coalescent events

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