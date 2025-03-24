using Test, PhyloNetworks
using InPhyNet

@testset "Main InPhyNet algo" begin

    ## Relatively simple test w/ a tree
    T = readnewick("((((A,B),C),(D,E)),(((F,G),(H,I)),((J,(K,L)),((M,N),O))));");
    constraints = [
        readnewick("(((A,B),C),(H,I));"),
        readnewick("((D,E),((M,N),O));"),
        readnewick("((F,G),(J,(K,L)));")
    ];
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
    ]);
    for i=1:size(D, 1) for j=(i+1):size(D,1) D[j,i] = D[i,j] end end
    mnet = inphynet(D, constraints, [l.name for l in T.leaf])
    @test hardwiredclusterdistance(T, mnet, false) == 0




    # Simple example w/ a network
    N = readnewick("(((A,(B)#H1),((C,#H1),D)),((E,(F)#H2),((G,#H2),H)));")
    constraints = [
        readnewick("((A,(B)#H1),((C,#H1),D));");
        readnewick("((E,(F)#H2),((G,#H2),H));")
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
    mnet = inphynet(D, constraints, sort([l.name for l in N.leaf]))
    @test hardwiredclusterdistance(N, mnet, false) == 0


    # More complicated example w/ a network
    N = readnewick("((((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8))),((t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17))))))),(((t18,t19),(t20,(t21,(t22,(t23,t24))))),((((t25,t26),(t27,t28)),(((t29,t30),(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,((t37,t38),(t39,t40)))))));")
    constraints = [     # N was built from the constraints as (c[1],(c[2],(c[3],c[4])))
        readnewick("(((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8)));"),
        readnewick("(t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17)))))));"),
        readnewick("((t18,t19),(t20,(t21,(t22,(t23,t24)))));"),
        readnewick("((((t25,t26),(t27,t28)),(((t29,t30),(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,((t37,t38),(t39,t40))));")
    ]
    D, names = majorinternodecount(N)
    mnet = inphynet(D, constraints, names)
    @test hardwiredclusterdistance(mnet, N, false) == 0


    # Same example as above but duplicated into (A,A)
    N = readnewick("(((((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8))),((t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17))))))),(((t18,t19),(t20,(t21,(t22,(t23,t24))))),((((t25,t26),(t27,t28)),(((t29,t30),(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,((t37,t38),(t39,t40))))))),((((t41,(t42)#H5),((t43,#H5),t44)),((t45,(t46)#H6),((t47,#H6),t48))),((t49,(t50,((t51,#H7),(t52,(((t53,t54))#H7,(t55,(t56,t57))))))),(((t58,t59),(t60,(t61,(t62,(t63,t64))))),((((t65,t66),(t67,t68)),(((t69,t70),(t71,t72)),#H8)),((((t73,t74),(t75,t76)))#H8,((t77,t78),(t79,t80))))))));")
    constraints = [     # N was built from the constraints as (c[1],(c[2],(c[3],c[4])))
        readnewick("(((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8)));"),
        readnewick("(t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17)))))));"),
        readnewick("((t18,t19),(t20,(t21,(t22,(t23,t24)))));"),
        readnewick("((((t25,t26),(t27,t28)),(((t29,t30),(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,((t37,t38),(t39,t40))));"),
        readnewick("(((t41,(t42)#H1),((t43,#H1),t44)),((t45,(t46)#H2),((t47,#H2),t48)));"),
        readnewick("(t49,(t50,((t51,#H3),(t52,(((t53,t54))#H3,(t55,(t56,t57)))))));"),
        readnewick("((t58,t59),(t60,(t61,(t62,(t63,t64)))));"),
        readnewick("((((t65,t66),(t67,t68)),(((t69,t70),(t71,t72)),#H4)),((((t73,t74),(t75,t76)))#H4,((t77,t78),(t79,t80))));")
    ]
    D, names = majorinternodecount(N)
    mnet = inphynet(D, constraints, names)
    @test hardwiredclusterdistance(mnet, N, false) == 0

end


@testset "InPhyNet - pruned subtrees with re-rooting" begin
    for m in [5, 10, 25, 50, 100]
        T = readnewick(joinpath(@__DIR__, "examples", "n250.tre"))
        T_subsets = centroid_edge_decomposition(T, 5, m)
        T_constraints = prune_network(T, T_subsets)
        for c in T_constraints
            PhyloNetworks.fuseedgesat!(c.rooti, c)
        end

        D, namelist = internodedistance(T)
        mnet = inphynet(D, T_constraints, namelist)
        @test hardwiredClusterDistance(T, mnet, false) == 0

        pair_mnet = InPhyNet.inphynet_pairwise(D, T_constraints, namelist)
        @test hardwiredClusterDistance(T, pair_mnet, false) == 0
    end
end


@testset "InPhyNet - n500 networks" begin
    nets = readmultinewick(joinpath(@__DIR__, "examples", "n500.net"))
    for net in nets
        PhyloNetworks.addleaf!(net, getroot(net), "OUTGROUP", 2.0)
        subsets = centroid_edge_decomposition(majortree(net), 5, 50)
        constraints = prune_network(net, subsets)
        D, namelist = majorinternodedistance(net)
        mnet = inphynet(D, constraints, namelist)
        rootatnode!(mnet, "OUTGROUP")
        @test hardwiredClusterDistance(mnet, net, true) <= 3
    end
end


@testset "InPhyNet - root reticulations" begin
    constraints = [
        readnewick("(#H1,(((((a)#H1,b),(c,d)),((e,f),(g,h)))));"),
        readnewick("(((i,j),k),(l,m));")
    ];
    try
        calculateAGID([majortree(c) for c in constraints])
        @test false
    catch e
        @test typeof(e) <: ArgumentError
    end
    init_mnet = inphynet(zeros(13, 13), constraints, vcat(tiplabels(constraints[1]), tiplabels(constraints[2])))

    D, namelist = majorinternodecount(bad_mnet)
    mnet = inphynet(D, constraints, namelist)

    @test hardwiredClusterDistance(mnet, init_mnet, false) == 0
end

