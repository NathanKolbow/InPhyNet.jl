# Good example net to work with: ((A,((B,(C,D)))#H1),(#H1,E));
# ---> _conditionOnCoalescences on the hybrid should give ? outputs (used to be 7, but now it's clear that's not the case)
# ---> 23 total parental trees (I believe, if output shows differently then may have to re-sketch)
include("./main.jl")

# ipt = IPT(readTopology("(10,(#H2,(1,(2,(((9)#H1,(3,(4,((5,6),(7,(8,#H1)))))))#H2))))root;"))
# _overwritemissinggammas!(top(ipt))

# ldict = LDict()
# i = 1
# for node in top(ipt).node
#     if node.leaf
#         ldict[node] = LineageNode(i)
#         i += 1
#     else
#         ldict[node] = nothing
#     end
# end

# output = _conditiononcoalescences(ipt, top(ipt).hybrid[1], ldict)
# # Looks good! (9 outputs)

# # output = _condensedivisions(output, [out.hybrid[1] for out in output], ldict)

# output = [_conditiononreticulation(out, top(out).hybrid[1], ldict) for out in output]
# for (j, set) in enumerate(output)
#     probsum = sum([prob(ipt) for ipt in set])
#     if probsum != 1.
#         println(set)
#         error("output["*string(j)*"] summed to "*string(probsum))
#     end
# end

# output = [_conditiononreticulation(tempnet, tempnet.hybrid[1], ldict) for tempnet in output]

ptrees, _ = getparentaltrees(readTopology("((A,((B,(C,D)))#H1),(#H1,E));"))
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees(readTopology("(((A,B),#H1), (((C,(D,#H2)))#H1,((E)#H2,F)));"))
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees("(10,(#H2,(1,(2,(((9)#H1,(3,(4,((5,6),(7,(8,#H1)))))))#H2))))root;")
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")
# Need coal probs to cover (N,O) := (4,1)



ipt = IPT(readTopology("(10,(#H2,(1,(2,(((9)#H1,(3,(4,((5,6),(7,(8,#H1)))))))#H2))))root;"))
_overwritemissinggammas!(top(ipt))
_overwritemissingbranchlengths!(top(ipt))

ldict = LDict()
i = 1
for node in top(ipt).node
    if node.leaf
        ldict[node] = LineageNode(i)
        i += 1
    else
        ldict[node] = nothing
    end
end

node, nodeidx = _getnexthybrid(top(ipt))
divisions = _conditiononcoalescences(ipt, node, ldict)
ipt = divisions[1]
node, nodeidx = _getnexthybrid(top(ipt))
net = top(ipt)
# After this point...

# Split reticulation seems to work as expected
splitreticulation!(net, node, findfirst(net.node .== [node]), LineageNode(Lineage(1)), LineageNode(Lineage(2)), ldict)
findall([node.name == "s_H1" for node in net.node])
ldict[net.node[5]]
ldict[net.node[24]]

# But conditioning on retic then coalescence breaks things
findall([node.name == "H1" for node in net.node])
divisions1 = _conditiononreticulation(ipt, node, ldict)         # KEY ERROR
ipt = divisions1[1]     # divisions1[2] *might* work, but runs into coalescent probability error [(N,O) := (4,1)]
net = top(ipt)
findall([node.name == "H1" for node in net.node])
findall([node.name == "s_H1" for node in net.node])     # only 1 result for [1] and [2] ????
node, nodeidx = _getnexthybrid(top(ipt))
divisions2 = _conditiononcoalescences(ipt, node, ldict)


# registered in `ldict` from the `net.node` list
listnode = net.node[5]
ldict[listnode]

# registered in `ldict` when following the topology
fromtop = getchildren(getchildren(getchildren(getchildren(getchildren(getchildren(getchildren(net.hybrid[1])[1])[1])[2])[2])[2])[2])[2]
ldict[fromtop]

# but these are NOT the same node, which is likely the source of the error
fromtop == listnode





badnode = top(ipt).node[5]                                      # only 1 node remains from the split; `s_H1`
badnode2 = getchildren(getchildren(getchildren(getchildren(getchildren(getchildren(getchildren(node)[1])[1])[2])[2])[2])[2])[2]
getchildren(badnode)
getchildren(badnode2)

println(repr(UInt64(pointer_from_objref(top(ipt).node[5]))))