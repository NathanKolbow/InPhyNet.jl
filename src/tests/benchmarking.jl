#############################################################################
# For benchmarking various code changes to see if speed improves or worsens #
#############################################################################
include("../main.jl")
using BenchmarkTools

#####################
# Current benchmark #
#####################

function benchfunc(complog::Bool)
    ipt, ldict, _ = _initparentaltreealgo(readTopology("(10,(#H2,(1,(2,(((9)#H1,(3,(6,(8,#H1)))))#H2))))root;"), neverwarn=true)
    _ = _getparentaltrees(ipt, ldict, usecomplog=complog)
end

@benchmark benchfunc(false)     # 238ms +- 16.ms (15.34 MiB, 196292 allocs)
@benchmark benchfunc(true)      # 229ms +- 3.3ms (14.76 MiB, 185945 allocs)


####################
# Previous results #
####################

# Using `complog` in `_calculatetotalcoalescentprobability`
#
# Results:
# - without:    238ms +- 16.ms (15.34 MiB, 196292 allocs)
# - with:       229ms +- 3.3ms (14.76 MiB, 185945 allocs)