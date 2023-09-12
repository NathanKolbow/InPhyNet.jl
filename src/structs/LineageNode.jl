# We use these nodes to hold the information that, at the point of the node,
# each of the lineages in `lineages` are present and *have not coalesced yet*
mutable struct LineageNode <: PhyloNetworks.ANode
    edge::Vector{PhyloNetworks.EdgeT{<:PhyloNetworks.ANode}}
    lineages::Vector{Lineage}
end

nlineages(node::LineageNode) = length(lineages)
lineages(node::LineageNode) = node.lineages