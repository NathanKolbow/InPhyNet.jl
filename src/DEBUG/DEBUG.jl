######################################################
# We trace a single standalone problem in this file. #
######################################################
include("./DEBUG-SETUP.jl")


lnode = LineageNode([Lineage(1), Lineage(2), Lineage(3), Lineage(4)])
opts, probs = getcoalescentcombos(lnode, 1.)
sum(probs)
length(opts)

1 + binomial(4,2) +
    binomial(4,2)*binomial(3,2) +
    binomial(4,2)*binomial(3,2)*1
