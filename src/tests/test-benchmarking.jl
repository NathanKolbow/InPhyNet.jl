#############################################################################
# For benchmarking various code changes to see if speed improves or worsens #
#############################################################################

include("../main.jl")
using BenchmarkTools

#####################
# Current benchmark #
#####################

function benchfunc()
    ipt, ldict, _ = _initparentaltreealgo(readTopology("(10,(#H2,(1,(2,(((9)#H1,(3,(6,(8,#H1)))))#H2))))root;"), neverwarn=true)
    _ = _getparentaltrees(ipt, ldict)
end

# Pre-computed binoms
@benchmark benchfunc()     # 143ms +- 4.9ms (8.74 MiB, 117499 allocs)
# Not pre-computed binoms
@benchmark benchfunc()     # 142ms +- 3.8ms (8.78 MiB, 117509 allocs)


####################
# Previous results #
####################

# Using `complog` in `_calculatetotalcoalescentprobability`
#
# Results:
# - without:    238ms +- 16.ms (15.34 MiB, 196292 allocs)
# - with:       229ms +- 3.3ms (14.76 MiB, 185945 allocs)