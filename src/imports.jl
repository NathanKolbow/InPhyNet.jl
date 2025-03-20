using PhyloNetworks, DataStructures, Graphs, StatsBase

import PhyloNetworks: deleteNode!, deleteEdge!, addhybridedge!, fuseedgesat!
# import Combinatorics: combinations, partitions, powerset
import Base: ==, *, -, names, getindex, setindex!, @warn
import Graphs: SimpleGraph, a_star, add_edge!