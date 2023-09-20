# Good example net to work with: ((A,((B,(C,D)))#H1),(#H1,E));
# ---> _conditionOnCoalescences on the hybrid should give ? outputs (used to be 7, but now it's clear that's not the case)
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

output = _conditiononcoalescences(net, net.hybrid[1], ldict)
# Looks good! (9 outputs)

output = _condensedivisions(output, [out.hybrid[1] for out in output], ldict)

_conditiononreticulation(output[1], output[1].hybrid[1], ldict)

output = [_conditiononreticulation(tempnet, tempnet.hybrid[1], ldict) for tempnet in output]

ptrees, ldict = getparentaltrees(readTopology("((A,((B,(C,D)))#H1),(#H1,E));"))
t = ptrees[1]   # all lineages go RIGHT in this tree
n = t.node[t.root]
ldict[getchildren(getchildren(n)[1])[2]]    # this USED to lead to the stale SPLIT-H1 node, but now that doesn't exist
ldict[getchildren(getchildren(n)[2])[1]]    # all lineages



# for testing splitreticulation:
# (((A,B),#H1), (((C,(D,#H2)))#H1,((E)#H2,F)));
# H2 is below H1