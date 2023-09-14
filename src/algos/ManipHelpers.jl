import DataStructures: Queue
# Helper functions for manipulating PhyloNetworks network structures

function copyldictcontents!(oldnet::HybridNetwork, newnet::HybridNetwork, ldict::LDict)
    for (newnode, oldnode) in zip(newnet.node, oldnet.node)
        lineagedict[newnode] = lineagedict[oldnode]
    end
end

# Splits the given reticulation into 2 tree-like edges
# Creates a deep copy of the net in the process
function splitreticulation(net::HybridNetwork, retic::PhyloNetworks.Node, leftline::LineageNode, rightline::LineageNode, ldict::LDict)

end