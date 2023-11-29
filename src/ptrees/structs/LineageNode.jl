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

const LDict = Dict{PhyloNetworks.Node, Union{LineageNode, Nothing}}

# Pretty-printing so that LineageNode is quickly and easily interpretable
Base.show(io::IO, x::LineageNode) = print(io, condense(x.lineages))
Base.show(io::IO, m::MIME"text/plain", x::LineageNode) = print(io, prettyformat([condense(l) for l in x.lineages]))

function prettyformat(outputs)
    ret = ""
    for (i, output) in enumerate(outputs)
        ret *= replace(string(output), "Any" => "", "Vector" => "", "{" => "", "}" => "")
        if i != length(outputs)
            ret *= "\n"
        end
    end
    return ret
end

condense(l1::Lineage, l2::Lineage) = condense([l1, l2])

function Base.:(==)(l1::LineageNode, l2::LineageNode)
    if length(lineages(l1)) != length(lineages(l2)) return false end
    if length(lineages(l1)) == 1 return lineages(l1)[1] == lineages(l2)[1] end

    l2idxs = collect(eachindex(lineages(l2)))
    for i in eachindex(lineages(l1))
        if i == length(lineages(l1)) break end
        foundequality = -1
        for j in l2idxs
            if lineages(l1)[i] == lineages(l2)[j]
                foundequality = j
                break
            end
        end
        if foundequality == -1 return false end
        deleteat!(l2idxs, foundequality)
    end

    return true
end


ln1 = LineageNode(Lineage(Lineage(Lineage(1), Lineage(2)), Lineage(3)))
ln2 = LineageNode(Lineage(Lineage(3), Lineage(Lineage(1), Lineage(2)))) 
ln3 = LineageNode(Lineage(Lineage(1), Lineage(Lineage(3), Lineage(2))))

ln1 == ln2 || error("LineageNode equality test failed.")
ln1 != ln3 || error("LineageNode equality test failed.")
ln2 != ln3 || error("LineageNode equality test failed.")