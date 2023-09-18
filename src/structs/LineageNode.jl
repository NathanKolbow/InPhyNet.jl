using PhyloNetworks
import Combinatorics: partitions

# We use these nodes to hold the information that, at the point of the node,
# each of the lineages in `lineages` are present and *have not coalesced yet*
mutable struct LineageNode <: PhyloNetworks.ANode
    # Fields
    lineages::Vector{Lineage}
    # children::Vector{LineageNode}

    # Constructors
    function LineageNode(l1::LineageNode, l2::LineageNode)
        new(vcat(lineages(l1), lineages(l2)))
    end
    function LineageNode(i::Integer)
        new(Vector{Lineage}([Lineage(i)]))
    end
    LineageNode(l::Lineage) = new([l])
    function LineageNode(ls::Vector{Lineage})
        new(ls)
    end
    function LineageNode()
        new([])
    end
end

# Helper methods
nlineages(node::LineageNode) = length(node.lineages)
nlineages(node::Nothing) = 0
lineages(node::LineageNode) = node.lineages

# Gets all of the possible coalescent combinations for
# the lineages in `l.lineages`
function getcoalescentcombos(l::LineageNode)
    parts = partitions(lineages(l))
    returnlist = Array{LineageNode}(undef, length(parts))

    for (i, sets) in enumerate(parts)
        templin = LineageNode()
        for set in sets
            coal = coalesce(set)
            if typeof(coal) <: LineageNode
                push!(lineages(templin), coalesce(set))
            else
                for co in coal push!(lineages(templin), co) end
            end
        end
        returnlist[i] = templin
    end

    return returnlist
end

const LDict = Dict{PhyloNetworks.Node, Union{LineageNode, Nothing}}

# Pretty-printing so that LineageNode is quickly and easily interpretable
Base.show(io::IO, x::LineageNode) = print(io, condense(x.lineages))
Base.show(io::IO, m::MIME"text/plain", x::LineageNode) = print(io, prettyformat([condense(l) for l in x.lineages]))

function prettyformat(outputs)
    ret = ""
    for (i, output) in enumerate(outputs)
        ret *= replace(replace(replace(replace(string(output), "Any" => ""), "Vector" => ""), "{" => ""), "}" => "")
        if i != length(outputs)
            ret *= "\n"
        end
    end
    return ret
end

function condense(ls)
    if typeof(ls) <: Lineage 
        return condense(ls.lineage)
    elseif typeof(ls) <: Vector && length(ls) == 1
        if typeof(ls[1]) <: Lineage
            return condense(ls[1])
        end
        return ls[1]
    end
    return [condense(l) for l in ls]
end

condense(l1::Lineage, l2::Lineage) = condense([l1, l2])