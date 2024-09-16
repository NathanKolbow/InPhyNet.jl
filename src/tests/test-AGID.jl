using InPhyNet, PhyloNetworks, NamedArrays, Test

gt1 = readTopology("((A:1,B:1):1,(C:1,D:1):1);")
D, namelist = internodedistance(gt1)
n = NamedArray(D, (namelist, namelist))

@test n["A", "B"] == 2
@test n["C", "D"] == 2
@test n["A", "C"] == 4
@test n["A", "D"] == 4
@test n["B", "C"] == 4
@test n["B", "D"] == 4


gt2 = readTopology("((A:2,B:2):1,(C:2,D:2):1);")
D, namelist = internodedistance(gt2)
n = NamedArray(D, (namelist, namelist))

@test n["A", "B"] == 4
@test n["C", "D"] == 4
@test n["A", "C"] == 6
@test n["A", "D"] == 6
@test n["B", "C"] == 6
@test n["B", "D"] == 6

