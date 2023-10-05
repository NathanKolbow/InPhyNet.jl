# Good example net to work with: ((A,((B,(C,D)))#H1),(#H1,E));
# ---> _conditionOnCoalescences on the hybrid should give ? outputs (used to be 7, but now it's clear that's not the case)
# ---> 23 total UNIQUE parental trees (fairly certain, not 100% though)
include("./main.jl")


ptrees, _ = getparentaltrees(readTopology("((A,((B,(C,D)))#H1),(#H1,E));"))
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees(readTopology("(((A,B),#H1), (((C,(D,#H2)))#H1,((E)#H2,F)));"))
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees("(10,(#H2,(1,(2,(((9)#H1,(7,(8,#H1))))#H2))))root;")
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees("(10,(#H2,(1,(2,(((9)#H1,(3,(6,(8,#H1)))))#H2))))root;")
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees("(10,(#H2,(1,(2,(((9)#H1,(3,(4,(5,((6,7),(8,#H1)))))))#H2))))root;")
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")
# Need coal probs to cover (N,O) := (5,1) and (6,1)



# Computation log efficiency improvement test
using BenchmarkTools

function benchfunc(complog::Bool)
    ipt, ldict, _ = _initparentaltreealgo(readTopology("(10,(#H2,(1,(2,(((9)#H1,(3,(6,(8,#H1)))))#H2))))root;"), neverwarn=true)
    _ = _getparentaltrees(ipt, ldict, usecomplog=complog)
end

@benchmark benchfunc(false)     # 238ms +- 16.ms (15.34 MiB, 196292 allocs)
@benchmark benchfunc(true)      # 229ms +- 3.3ms (14.76 MiB, 185945 allocs)