# Good example net to work with: ((A,((B,(C,D)))#H1),(#H1,E));
# ---> _conditionOnCoalescences on the hybrid should give 7 output
# ---> 23 total parental trees (I believe, if output shows differently then may have to re-sketch)
include("./main.jl")

net = readTopology("((A,((B,(C,D)))#H1),(#H1,E));")

ldict = Dict{PhyloNetworks.Node, Union{LineageNode, Nothing}}()
i = 1
for node in net.node
    if node.leaf
        ldict[node] = LineageNode(i, net.numTaxa)
        i += 1
    else
        ldict[node] = nothing
    end
end

output = _conditionOnCoalescences(net, net.hybrid[1], ldict)
# Looks good!

ptrees = getParentalTrees(readTopology("((A,((B,(C,D)))#H1),(#H1,E));"))