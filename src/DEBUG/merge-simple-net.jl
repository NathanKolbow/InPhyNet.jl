include("../main.jl")
# debugging simple network example
N = readTopology("(((A,B)#H1,((C,#H1),D)),((E,F)#H2,((G,#H2),H)));")
constraints = [
    readTopology("((A,B)#H1,((C,#H1),D));");
    readTopology("((E,F)#H2,((G,#H2),H));")
]
D = Matrix{Float64}([   # major tree internode distances
    0.  1.  3.  3.  5.  5.  5.  5.;
    -1. 0.  3.  3.  5.  5.  5.  5.;
    -1. -1. 0.  1.  5.  5.  5.  5.;
    -1. -1. -1. 0.  5.  5.  5.  5.;
    -1. -1. -1. -1. 0.  1.  3.  3.;
    -1. -1. -1. -1. -1. 0.  3.  3.;
    -1. -1. -1. -1. -1. -1. 0.  1.;
    -1. -1. -1. -1. -1. -1. -1. 0.;
])
for i=1:size(D, 1) for j=(i+1):size(D,1) D[j,i] = D[i,j] end end

# now this is the inside of njmerge
n = size(D, 1)
names = sort([l.name for l in N.leaf])
subnets = Vector{SubNet}([SubNet(i, names[i]) for i in 1:n])
reticmap = ReticMap(constraints)

# infinite loop on second run
##########################################################################################################################################################
# while n > 2
possible_siblings = findvalidpairs(constraints, names)

# Find optimal (i, j) idx pair for matrix Q
i, j = findoptQidx(D, possible_siblings)

# connect subnets i and j
# TODO: before implementing this section, sketch out
#       what it should look like to connect 2 subnets
#       (review nj! code)
subnets[i], edgei, edgej = mergesubnets!(subnets[i], subnets[j])
updateconstraints!(names[i], names[j], constraints, reticmap, edgei, edgej)

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
names = view(names, idxfilter)

n -= 1
##########################################################################################################################################################
mnet = mergesubnets!(subnets[1], subnets[2])[1]
mnet = HybridNetwork(mnet)
placeretics!(mnet, reticmap)