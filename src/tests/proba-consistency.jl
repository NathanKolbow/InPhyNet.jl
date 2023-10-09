include("../main.jl")

"""
Tests the above probabilities for `N` input lineages
and branch length `bl`.
"""
function probatest(N::Int64, bl::Real)
    probsum = 0.
    for nout=1:N
        prob = _calculatecoalescentprobability(N, nout, bl)
        mult = prod([binomial(i, 2) for i=(nout+1):N])
        probsum += mult * prob
    end
    abs(probsum - 1.) < 1e-12 || error("probsum: "*string(probsum))
end

probatest(4, 1.) || error("failed")
probatest(5, 1.) || error("failed")
probatest(6, 1.) || error("failed")
probatest(7, 1.) || error("failed")
probatest(8, 1.) || error("failed")

"""
Tests the above probabilities for `N` input lineages
and branch length `bl`.
"""
function probatoptest(N::Int64, bl::Real)
    N > 0 || error("N>0 required.")
    bl = BigFloat(bl)
    lnode = LineageNode([Lineage(i) for i=1:N])
    combos, probs = getcoalescentcombos(lnode, bl)

    println(combos)
    println([round(Float64(p), digits=5) for p in probs])
    abs(sum(probs) - BigFloat(1.)) < 1e-12 || error("sum(probs) = "*string(sum(probs)))
end

probatoptest(3, 1.)
probatoptest(4, 1.)
probatoptest(5, 1.)
probatoptest(6, 1.)
probatoptest(7, 1.)