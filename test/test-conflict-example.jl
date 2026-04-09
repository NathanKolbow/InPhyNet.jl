loadconflictexample()

@test_throws "No compatible merge found." inphynet(D, constraints, namelist; refuse_pairwise=true)
@test isa(inphynet(D, constraints, namelist), HybridNetwork)
@test isa(InPhyNet.inphynetpairwise(D, constraints, namelist), HybridNetwork)