using PhyloNetworks, DataStructures, DataFrames, CSV, Graphs

import PhyloNetworks: deleteNode!, deleteEdge!, addhybridedge!, fuseedgesat!
import Combinatorics: combinations
import Base: ==, *, -, names, getindex, setindex!
import Combinatorics: partitions, powerset
import Test: @warn
import Graphs: SimpleGraph, a_star, add_edge!


# Just for dev, to be removed later
using Plots, GraphRecipes, PhyloCoalSimulations, PhyloPlots