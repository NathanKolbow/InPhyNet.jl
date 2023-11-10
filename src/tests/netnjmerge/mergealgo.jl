include("../../main.jl")

# Test w/ a tree example
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
mergednet = netnj!(D, constraints, names=T.names)
mergednet = HybridNetwork(mergednet.nodes, mergednet.edges)
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
end
mindist
mindist == 0 || error("")