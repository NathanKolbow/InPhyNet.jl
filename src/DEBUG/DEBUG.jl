using Revise
include("../../simulations/pipelines.jl")
include("../../simulations/robustness-fxns.jl")




# findvalidpairs issue
cs = [
    readTopology("(t1,t5);"),
    readTopology("(t1,(t5,t7));")
]
findvalidpairs(cs, cs[2].names)
findvalidpairs([cs[2]], cs[2].names)




####
truenet, cs = loadTrueData("n40h4")
D, namelist = majorinternodedistance(truenet)
D .= 0.
mnet = netnj!(D, cs, namelist)

namelistNew = ["t1", "t11", "t13", "t16"]
possible_siblings = NetMerge.findvalidpairs(constraints, namelistNew)




# Looks like constraints[2] is the one that's causing the problem
# - topology updates are bad at some point; constraints[2] has 4 tips but only 1 shows
c2 = readTopology("(t9,(t10,((t11,#H3)intA,(t12,(((t13,t14)intB)#H3,(t15,(t16,t17)intC)intD)intE)intF)intG)intH)r;")
# PhyloPlots.plot(c2)

r = ReticMap([c2])
e() = Edge(-1)

updateconstraints!("t1", "t2", [c2], r, e(), e())
updateconstraints!("t1", "t9", [c2], r, e(), e())
updateconstraints!("t1", "t18", [c2], r, e(), e())
updateconstraints!("t1", "t19", [c2], r, e(), e())
updateconstraints!("t1", "t25", [c2], r, e(), e())
updateconstraints!("t1", "t26", [c2], r, e(), e())
updateconstraints!("t3", "t4", [c2], r, e(), e())
updateconstraints!("t1", "t3", [c2], r, e(), e())
updateconstraints!("t5", "t6", [c2], r, e(), e())
updateconstraints!("t5", "t10", [c2], r, e(), e())
updateconstraints!("t5", "t20", [c2], r, e(), e())
updateconstraints!("t5", "t27", [c2], r, e(), e())
updateconstraints!("t5", "t28", [c2], r, e(), e())
updateconstraints!("t7", "t8", [c2], r, e(), e())
updateconstraints!("t5", "t7", [c2], r, e(), e())
updateconstraints!("t1", "t5", [c2], r, e(), e())   # issue emerges here
c2





updateconstraints!("t11", "t12", [c2], r, e(), e())
updateconstraints!("t13", "t14", [c2], r, e(), e())
updateconstraints!("t13", "t15", [c2], r, e(), e())
updateconstraints!("t16", "t17", [c2], r, e(), e())

# minimal example
net = readTopology("((E,((A,#H1),(((B)#H1,C),D))),(F,(G,(H,I))));")
r = ReticMap([net])
updateconstraints!("A", "B", [net], r, e(), e())
net

# constraints[4] is also a problem
# - shows negative edges
