import Combinatorics: partitions

# We use these nodes to hold the information that, at the point of the node,
# each of the lineages in `lineages` are present and *have not coalesced yet*
mutable struct LineageNode <: PhyloNetworks.ANode
    # Fields
    lineages::Vector{Lineage}
    children::Vector{LineageNode}

    # Constructors
    function LineageNode(l1::LineageNode, l2::LineageNode)
        new(vcat(lineages(l1), lineages(l2)))
    end
    function LineageNode(i::Integer, ntaxa::Integer)
        new(Vector{Lineage}([Lineage(i, ntaxa)]))
    end
    function LineageNode(ls::Vector{Lineage})
        new(ls)
    end
    function LineageNode()
        new([])
    end
end

# Helper methods
nlineages(node::LineageNode) = length(lineages)
lineages(node::LineageNode) = node.lineages

# Gets all of the possible coalescent combinations for
# the lineages in `l.lineages`
function getcoalescentcombos(l::LineageNode)
    parts = partitions(lineages(l))
    returnlist = Array{LineageNode}(undef, length(parts))

    for (i, sets) in enumerate(parts)
        templin = LineageNode()
        for set in sets
            push!(lineages(templin), coalesce(set))
        end
        returnlist[i] = templin
    end

    return returnlist
end