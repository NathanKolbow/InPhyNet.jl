include("../../main.jl")

# Complicated example w/ a network
truenet = readTopology("((((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8))),((t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17))))))),((t18,(t20,(t21,(t22,t23)))),(((t25,t27),((t29,(t31,t32)),#H4)),((((t33,t34),(t35,t36)))#H4,(t37,t38))))));")
ngt = 1000

# 1. simulate gene trees & save them
for e in truenet.edge
    e.length = 1.
    if e.hybrid
        e.gamma = 0.5
    end
end
sims = simulatecoalescent(truenet, ngt, 1)
simfile = tempname()
writeMultiTopology(sims, simfile)

# 2. read the sims and run MSCquartets
ptablefile = tempname()
ptablefile = "src/tests/netnjmerge/1000currsim.hbptab"
# in mscquartets.R:
#   write.HBpTable($simfile, $ptablefile)
# RCall is horrendously slow...

# 3. parse MSCquartets output & perform subset decomposition
hybsubsets, treesubset = decomposeFromQuartets(ptablefile, cutoff=0.000000000000000000000000000001)
# rm(simfile)
# rm(ptablefile)

# 4. pretend we got exact network estimates for each hybsubset and treesubset
constraints = [
    readTopology("(((t25,t27),((t31,t32),#H4)),((((t33,t34),(t35,t36)))#H4,(t37,t38)));"),
    readTopology("((t11,#H3),((t12,(t13,t14)))#H3);"),
    readTopology("((t1,(t2)#H1),((t3,#H1),t4));"),
    readTopology("((t5,(t6)#H2),((t7,#H2),t8));"),
    readTopology("((((t9,t10),(t15,(t16,t17))),(t18,(t20,(t21,(t22,t23))))),t29);")
]

# 5. AGID from true gene trees
D, namelist = calculateAGID(sims)

# 6. run the algo!
mnet = netnj!(D, constraints, namelist=namelist)

hardwiredClusterDistance(mnet, truenet, false)
# 8!