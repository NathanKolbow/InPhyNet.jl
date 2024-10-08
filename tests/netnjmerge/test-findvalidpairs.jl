# All necessary portions of the package should already be imported
constraints = [
    readTopology("((A,((B)#H2)#H1), (((C,#H2),#H1),D));"),
    readTopology("(((A)#H1,(E,#H1)),((F)#H2,(G,#H2)));"),
    readTopology("((B,E),(G,H));"),
    readTopology("((A,F),(G,I));")
]

leafnames = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]

# constraints[1] only
out = findvalidpairs(constraints[1], leafnames)
all(out[5:10] .== 1) || error("test failed")
outview = view(out, 1:4, 1:4)
for (i, j) in [(1, 2), (2, 3), (3, 4), (2, 1), (3, 2), (4, 3)]
    outview[i, j] -= 1
end
all(outview .== 0) || error("test failed")

# constraints[2] only
out = findvalidpairs(constraints[2], leafnames[[1, 5, 6, 7]])
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 3)]
    out[i, j] -= 1
end
all(out .== 0) || error("test failed")

# constraints[3] only
out = findvalidpairs(constraints[3], leafnames[[2, 5, 7, 8]])
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 3)]
    out[i, j] -= 1
end
all(out .== 0) || error("test failed")

# constraints[4] only
out = findvalidpairs(constraints[4], leafnames[[1, 6, 7, 9]])
for (i, j) in [(1, 2), (2, 1), (3, 4), (4, 3)]
    out[i, j] -= 1
end
all(out .== 0) || error("test failed")

# all constraints
out = findvalidpairs(constraints, leafnames)
for (i, j) in [(1, 2), (2, 3), (3, 4), (1, 5), (6, 7), (2, 5), (7, 8), (7, 9),
    (1, 8), (2, 6), (2, 9), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (5, 9), (6, 8), (9, 8)]
    out[i, j] -= 1
    out[j, i] -= 1
end
out[4,5:9] .-= 1
out[5:9,4] .-= 1
all(out[1:9, 1:9] .== 0) || error("test failed")
all(out[1:10,10] .== 1) && all(out[10,1:10] .== 1) || error("test failed")