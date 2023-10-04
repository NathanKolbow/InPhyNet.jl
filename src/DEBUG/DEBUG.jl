######################################################
# We trace a single standalone problem in this file. #
# ptrees, _ = getparentaltrees("(10,(#H2,(1,(2,(((9)#H1,(3,(5,(7,(8,#H1))))))#H2))))root;")
######################################################
include("./DEBUG-SETUP.jl")

net = top(ipt)
for node in net.node
    ldict[node]     # everything in the input `net.node` is in the ldict
end

nodeQ = Queue{Node}()
enqueue!(nodeQ, net.node[net.root])

badnodes = []
while !isempty(nodeQ)
    currnode = dequeue!(nodeQ)
    try
        ldict[currnode]
        if !(currnode in net.node) push!(badnodes, currnode) end
    catch e
        push!(badnodes, currnode)
    end

    if !currnode.leaf 
        for child in getchildren(currnode) enqueue!(nodeQ, child) end
    end
end
length(badnodes)
badnodes[1] == badnodes[2]  # same node appeared twice b/c there is still a hybrid in the network
badnodes[3] == badnodes[4]  # same node appeared twice again b/c of the hybrid
badnodes[1] == badnodes[3]  # these nodes are different ????
# this is NOT the node causing the error, though...


nodeQ = Queue{Node}()
enqueue!(nodeQ, node)   # `node` is the next hybrid node being evaluated

badnodes = []
while !isempty(nodeQ)
    currnode = dequeue!(nodeQ)
    try
        ldict[currnode]
        if !(currnode in net.node) push!(badnodes, currnode) end
    catch e
        push!(badnodes, currnode)
    end

    if !currnode.leaf 
        for child in getchildren(currnode) enqueue!(nodeQ, child) end
    end
end
length(badnodes)    # we found s_H1 this time! it's not in `net.node`
badnodes[1]     # found by following topology, not in `net.node`
badnodes[2]     # in `net.node`, not in `ldict`

# KeyError
divisions = _conditiononcoalescences(ipt, node, ldict)




# Minimal example
include("../main.jl")
ipt, ldict, _ = _initparentaltreealgo(readTopology("(10,(#H2,(1,(2,(((9)#H1,(3,(5,(7,(8,#H1))))))#H2))))root;"))
ipt = _conditiononcoalescences(ipt, _getnexthybrid(top(ipt))[1], ldict)[1]
ipt = _conditiononreticulation(ipt, _getnexthybrid(top(ipt))[1], ldict)[1]

net = top(ipt)
nodeQ = Queue{Node}()
enqueue!(nodeQ, net.node[net.root])

badnodes = []
while !isempty(nodeQ)
    currnode = dequeue!(nodeQ)
    try
        ldict[currnode]
        if !(currnode in net.node) push!(badnodes, currnode) end
    catch e
        push!(badnodes, currnode)
    end

    if !currnode.leaf 
        for child in getchildren(currnode) enqueue!(nodeQ, child) end
    end
end
length(badnodes)




# Breaking down the example
ipt, ldict, _ = _initparentaltreealgo(readTopology("(10,(#H2,(1,(2,(((9)#H1,(3,(5,(7,(8,#H1))))))#H2))))root;"))
ipt = _conditiononcoalescences(ipt, _getnexthybrid(top(ipt))[1], ldict)[1]

net = top(ipt)
hyb = _getnexthybrid(top(ipt))[1]
child = ldict[hyb]
hybidx = findfirst(top(ipt).node .== [hyb])

newnet = deepcopy(net)
newhyb = newnet.node[hybidx]
copyldictcontents!(net, newnet, ldict)
majorline = LineageNode(lineages(child)[setdiff(1, 1)])     # majorline = [1], minorline = [setdiff] does NOT cause the same error.
minorline = LineageNode(lineages(child)[1])                     
splitreticulation!(newnet, newhyb, hybidx, majorline, minorline, ldict)

# Issues immediately following `splitreticulation!`
net = newnet
nodeQ = Queue{Node}()
enqueue!(nodeQ, net.node[net.root])   # `node` is the next hybrid node being evaluated

badnodes = []
while !isempty(nodeQ)
    currnode = dequeue!(nodeQ)
    try
        ldict[currnode]
        if !(currnode in net.node) push!(badnodes, currnode) end
    catch e
        push!(badnodes, currnode)
    end

    if !currnode.leaf 
        for child in getchildren(currnode) enqueue!(nodeQ, child) end
    end
end
length(badnodes)