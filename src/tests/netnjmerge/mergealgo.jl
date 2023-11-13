include("../../main.jl")

## Relatively simple test w/ a tree
T = readTopology("((((A,B),C),(D,E)),(((F,G),(H,I)),((J,(K,L)),((M,N),O))));")
constraints = [
    readTopology("(((A,B),C),(H,I));"),
    readTopology("((D,E),((M,N),O));"),
    readTopology("((F,G),(J,(K,L)));")
]
D = Matrix{Float64}([
     0.  1. 2.  4.  4.  8.  8.  8.  8.  8.  9.  9.  9.  9.  8.;
    -1.  0. 2.  4.  4.  8.  8.  8.  8.  8.  9.  9.  9.  9.  8.;
    -1. -1. 0.  3.  3.  6.  6.  6.  6.  7.  8.  8.  8.  8.  7.;
    -1. -1. -1. 0.  1.  7.  7.  7.  7.  8.  9.  9.  9.  9.  8.;
    -1. -1. -1. -1. 0.  7.  7.  7.  7.  8.  9.  9.  9.  9.  8.;
    -1. -1. -1. -1. -1. 0.  1.  3.  3.  6.  7.  7.  7.  7.  6.;
    -1. -1. -1. -1. -1. -1. 0.  3.  3.  6.  7.  7.  7.  7.  6.;
    -1. -1. -1. -1. -1. -1. -1. 0.  1.  6.  7.  7.  7.  7.  6.;
    -1. -1. -1. -1. -1. -1. -1. -1. 0.  6.  7.  7.  7.  7.  6.;
    -1. -1. -1. -1. -1. -1. -1. -1. -1. 0.  2.  2.  4.  4.  3.;
    -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. 0.  1.  5.  5.  4.;
    -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. 0.  5.  5.  4.;
    -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. 0.  1.  2.;
    -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. 0.  2.;
    -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. -1. 0.;
])
for i=1:size(D, 1) for j=(i+1):size(D,1) D[j,i] = D[i,j] end end
mergednet = netnj!(D, constraints, names=[l.name for l in T.leaf])
mergednet.root = mergednet.numNodes
mergednet = readTopology(writeTopology(mergednet))

writeTopology(mergednet)                        # identical up to root differences!
hardwiredClusterDistance(T, mergednet, true)    # not 0 b/c of root differences, silly...

mindist = Inf
minidx = -1
for i=1:mergednet.numNodes
    mergednet.root = i
    if hardwiredClusterDistance(T, mergednet, false) < mindist
        mindist = hardwiredClusterDistance(T, mergednet, false)
        minidx = i
    end
    mergednet.root = minidx
end
mindist
mindist == 0 || error("")


# Simple example w/ a network
N = readTopology("(((A,(B)#H1),((C,#H1),D)),((E,(F)#H2),((G,#H2),H)));")
constraints = [
    readTopology("((A,(B)#H1),((C,#H1),D));");
    readTopology("((E,(F)#H2),((G,#H2),H));")
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
mnet = netnj!(D, constraints, names=sort([l.name for l in N.leaf]))