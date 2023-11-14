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


# More complicated example w/ a network
N = readTopology("((((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8))),((t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17))))))),(((t18,t19),(t20,(t21,(t22,(t23,t24))))),((((t25,t26),(t27,t28)),(((t29,t30),(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,((t37,t38),(t39,t40)))))));")
constraints = [     # N was built from the constraints as (c[1],(c[2],(c[3],c[4])))
    readTopology("(((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8)));"),
    readTopology("(t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17)))))));"),
    readTopology("((t18,t19),(t20,(t21,(t22,(t23,t24)))));"),
    readTopology("((((t25,t26),(t27,t28)),(((t29,t30),(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,((t37,t38),(t39,t40))));")
]
D, names = majorinternodedistance(N)
mnet = netnj!(D, constraints, names=names)
hardwiredClusterDistance(mnet, N, false) == 0 || error("test failed")

# Same example as above but duplicated into (A,A)
N = readTopology("(((((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8))),((t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17))))))),(((t18,t19),(t20,(t21,(t22,(t23,t24))))),((((t25,t26),(t27,t28)),(((t29,t30),(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,((t37,t38),(t39,t40))))))),((((t41,(t42)#H5),((t43,#H5),t44)),((t45,(t46)#H6),((t47,#H6),t48))),((t49,(t50,((t51,#H7),(t52,(((t53,t54))#H7,(t55,(t56,t57))))))),(((t58,t59),(t60,(t61,(t62,(t63,t64))))),((((t65,t66),(t67,t68)),(((t69,t70),(t71,t72)),#H8)),((((t73,t74),(t75,t76)))#H8,((t77,t78),(t79,t80))))))));")
constraints = [     # N was built from the constraints as (c[1],(c[2],(c[3],c[4])))
    readTopology("(((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8)));"),
    readTopology("(t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17)))))));"),
    readTopology("((t18,t19),(t20,(t21,(t22,(t23,t24)))));"),
    readTopology("((((t25,t26),(t27,t28)),(((t29,t30),(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,((t37,t38),(t39,t40))));"),
    readTopology("(((t41,(t42)#H1),((t43,#H1),t44)),((t45,(t46)#H2),((t47,#H2),t48)));"),
    readTopology("(t49,(t50,((t51,#H3),(t52,(((t53,t54))#H3,(t55,(t56,t57)))))));"),
    readTopology("((t58,t59),(t60,(t61,(t62,(t63,t64)))));"),
    readTopology("((((t65,t66),(t67,t68)),(((t69,t70),(t71,t72)),#H4)),((((t73,t74),(t75,t76)))#H4,((t77,t78),(t79,t80))));")
]
D, names = majorinternodedistance(N)
mnet = netnj!(D, constraints, names=names)
hardwiredClusterDistance(mnet, N, false) == 0 || error("test failed")
