include("../main.jl")

# we want this constraint case to become ((A,B),C) OR (A,(B,C))
# while logging the reticulation to retic map
c = readTopology("((A,(B)#H1),((C,#H1),D));")
r = ReticMap([c])
e = Edge(1, -1.)

mergeconstraintnodes!(c, c.leaf[2], c.leaf[3], r, e, e)
c

# we want to retain this constraint case becoming (A,((C,#H1),D))
c = readTopology("((A,(B)#H1),((C,#H1),D));")
r = ReticMap([c])
e = Edge(1, -1.)

mergeconstraintnodes!(c, c.leaf[1], c.leaf[2], r, e, e)
c

# we want to retain this constraint case becoming (A,(C,D))
c = readTopology("(((A,#H1),(B)#H1),(C,D));")
r = ReticMap([c])
e = Edge(1, -1.)

mergeconstraintnodes!(c, c.leaf[1], c.leaf[2], r, e, e)
c


# retic map logging issue here
c = [readTopology("(((t25,t27),((t31,t32),#H4)),((((t33,t34),(t35,t36)))#H4,(t37,t38)));")]
D, namelist = majorinternodedistance(c[1])
n = size(D, 1)
r = ReticMap(c)
subnets = Vector{SubNet}([SubNet(i, namelist[i]) for i in 1:n])

while n > 6
    possible_siblings = findvalidpairs(c, namelist)
    
    # Find optimal (i, j) idx pair for matrix Q
    i, j = findoptQidx(D, possible_siblings)

    # connect subnets i and j
    subnets[i], edgei, edgej = mergesubnets!(subnets[i], subnets[j])
    updateconstraints!(namelist[i], namelist[j], c, r, edgei, edgej)

    # collapse taxa i into j
    for l in 1:n
        if l != i && l != j
            D[l, i] = D[i, l] = (D[l, i] + D[j, l] - D[i, j]) / 2
        end
    end

    # Remove data elements that corresponded to `j`
    idxfilter = [1:(j-1); (j+1):n]
    D = view(D, idxfilter, idxfilter)   # remove j from D
    subnets = view(subnets, idxfilter)
    namelist = view(namelist, idxfilter)

    n -= 1
end

# updateconstraints!("t25", "t37", c, r, edgei, edgej)
# here leads to a leaf w/ 0 edges
# follows path "c"
valpairs = findvalidpairs(c, namelist)
valpairs[1,6]

nodei = c[1].node[findfirst([n.name == "t25" for n in c[1].node])]
nodej = c[1].node[findfirst([n.name == "t37" for n in c[1].node])]

using Plots, GraphRecipes
graph = Graph(c[1], includeminoredges=true)
plot(graph)     # why is t25 --> t37

idxnodei = findfirst(c[1].node .== [nodei])
idxnodej = findfirst(c[1].node .== [nodej])
edgepath = a_star(graph, idxnodei, idxnodej)

nodesinpath = Array{Node}(undef, length(edgepath)+1)
edgesinpath = Array{Edge}(undef, length(edgepath))
for (i, gedge) in enumerate(edgepath)
    srcnode = c[1].node[gedge.src]
    dstnode = c[1].node[gedge.dst]

    if i == 1 nodesinpath[1] = c[1].node[gedge.src] end
    nodesinpath[i+1] = dstnode

    netedge = filter(e -> (srcnode in e.node) && (dstnode in e.node), dstnode.edge)[1]
    edgesinpath[i] = netedge
end

# find the node that should be the new tip after merging
newtip = nothing
for node in nodesinpath
    if node.leaf continue end
    isnewtip = true
    for edge in node.edge
        if edge.hybrid && !edge.isMajor
            isnewtip = false
            break
        end
    end
    if isnewtip
        newtip = node
        break
    end
end



# retic not getting logged
c = [readTopology("(((t25,t27),((t31,t32),#H4)),((((t33,t34),(t35,t36)))#H4,(t37,t38)));")]
namelist = [leaf.name for leaf in c[1].leaf]
n = c[1].numTaxa
subnets = Vector{SubNet}([SubNet(i, namelist[i]) for i in 1:n])
r = ReticMap(c)

# first round
mergeorder = [
    (3, 4),
    (7, 6),
    (5, 4),
    (5, 4),
    (1, 2),
    (4, 5),
    (1, 2),     # not logging properly here
    (3, 2),
    (2, 1)
]
for (i, j) in mergeorder
    println(c[1])
    
    subnets[i], edgei, edgej = mergesubnets!(subnets[i], subnets[j])
    updateconstraints!(namelist[i], namelist[j], c, r, edgei, edgej)
    idxfilter = [1:(j-1); (j+1):n]
    subnets = view(subnets, idxfilter)
    namelist = view(namelist, idxfilter)
    n -= 1
end
r