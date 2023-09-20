using PhyloNetworks
using DataStructures
import Base: ==
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
    partq = Queue{Any}()
    for part in parts enqueue!(partq, part) end

    returnlist = Vector{LineageNode}()
    while !isempty(partq)
        sets = dequeue!(partq)

        continuewhile = false
        for (i, set) in enumerate(sets)
            if length(set) > 2
                # If this grouping has > 2 taxa then we need to break it up further
                outcomes = coalesce(set)
                for (j, outcome) in enumerate(outcomes)
                    dupsets = ifelse(j == length(outcomes), sets, deepcopy(sets))
                    dupsets[i] = [outcome]
                    # re-queue the new copies
                    enqueue!(partq, dupsets)
                end

                # We've re-queued the updated objects, so move on in the queue
                continuewhile = true
                break
            end
        end

        if continuewhile continue end

        # If we didn't get caught in the above loop, that means this partition is good to be pushed
        templin = LineageNode()
        for lin in sets
            push!(lineages(templin), Lineage(lin)) 
        end
        push!(returnlist, templin)
    end

    # for (i, sets) in enumerate(parts)
    #     templin = LineageNode()

    #     for set in sets
    #         coal = coalesce(set)
    #         if typeof(coal) <: LineageNode
    #             push!(lineages(templin), coalesce(set))
    #         else
    #             for co in coal push!(lineages(templin), co) end
    #         end
    #     end
    #     returnlist[i] = templin
    # end

    return returnlist
end

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