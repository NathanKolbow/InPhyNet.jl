include("../main.jl")

# we want this constraint case to become ((A,B),C) OR (A,(B,C))
# while logging the reticulation to retic map
c = readTopology("((A,(B)#H1),((C,#H1),D));")
r = ReticMap([c])
e = Edge(1, -1.)

mergeconstraintnodes!(c, c.leaf[2], c.leaf[3], r, e, e)
c

# we want to retain this constraint case becoming (A,((C,#H1),D))
c = readTopology("((A,(B)#H1),((C,#H1),D));")
r = ReticMap([c])
e = Edge(1, -1.)

mergeconstraintnodes!(c, c.leaf[1], c.leaf[2], r, e, e)
c

# we want to retain this constraint case becoming (A,(C,D))
c = readTopology("(((A,#H1),(B)#H1),(C,D));")
r = ReticMap([c])
e = Edge(1, -1.)

mergeconstraintnodes!(c, c.leaf[1], c.leaf[2], r, e, e)
c