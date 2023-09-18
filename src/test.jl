# Good example net to work with: ((A,((B,(C,D)))#H1),(#H1,E));
# ---> _conditionOnCoalescences on the hybrid should give 7 output
# ---> 23 total parental trees (I believe, if output shows differently then may have to re-sketch)
include("./main.jl")

net = readTopology("((A,((B,(C,D)))#H1),(#H1,E));")

ldict = LDict()
i = 1
for node in net.node
    if node.leaf
        ldict[node] = LineageNode(i)
        i += 1
    else
        ldict[node] = nothing
    end
end

output = _conditionOnCoalescences(net, net.hybrid[1], ldict)
# Looks good! (should be 7)

output = [_conditionOnReticulation(tempnet, tempnet.hybrid[1], ldict) for tempnet in output]






ptrees = getParentalTrees(readTopology("((A,((B,(C,D)))#H1),(#H1,E));"))




# for testing splitreticulation:
# (((A,B),#H1), (((C,(D,#H2)))#H1,((E)#H2,F)));
# H2 is below H1