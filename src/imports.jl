using PhyloNetworks
using DataStructures
import Combinatorics: combinations
import Base: ==, *, -
import Combinatorics: partitions, powerset
import Test: @warn
import Graphs: SimpleGraph, a_star, add_edge!