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

@benchmark benchfunc()


####################
# Previous results #
####################

# Using `_combinedivisions!` after `_conditiononcoalescences`
# Results:
# - with:       98ms  +- 4.1ms (6.00 MiB, 82352 allocs)
# - without:    145ms +- 6.7ms (8.71 MiB, 117113 allocs)
#
# Using `complog` in `_calculatetotalcoalescentprobability`
#
# Results:
# - with:       229ms +- 3.3ms (14.76 MiB, 185945 allocs)
# - without:    238ms +- 16.ms (15.34 MiB, 196292 allocs)