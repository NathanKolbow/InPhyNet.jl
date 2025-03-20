using Test, PhyloNetworks
using InPhyNet
import InPhyNet: findvalidpairs, findsiblingpairs


@testset "Find valid pairs" begin

    # All necessary portions of the package should already be imported
    constraints = majortree.([
        readnewick("((A,(B)#H1), ((C,#H1),D));"),
        readnewick("(((A)#H1,(E,#H1)),((F)#H2,(G,#H2)));"),
        readnewick("((B,E),(G,H));"),
        readnewick("((A,F),(G,I));")
    ])


    leafnames = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]

    # constraints[1] only
    sibling_pairs = findsiblingpairs(constraints[1])
    out = findvalidpairs([constraints[1]], [sibling_pairs], leafnames)
    @test all(out[5:10] .== 1)
    outview = view(out, 1:4, 1:4)
    for (i, j) in [(1, 2), (3, 4), (2, 1), (4, 3)]
        outview[i, j] -= 1
    end
    @test all(outview .== 0)

    # all constraints
    sibling_pairs = [findsiblingpairs(c) for c in constraints]
    out = findvalidpairs(constraints, sibling_pairs, leafnames)
    for (i, j) in [(1, 2), (3, 4), (1, 5), (2, 5), (7, 8), (7, 9), (1, 8), (2, 6),
        (2, 9), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (5, 9), (6, 8), (9, 8)]
        out[i, j] -= 1
        out[j, i] -= 1
    end
    out[4,5:9] .-= 1
    out[5:9,4] .-= 1
    @test all(out[1:9, 1:9] .== 0)
    @test all(out[1:10,10] .== 1) && all(out[10,1:10] .== 1)

end

