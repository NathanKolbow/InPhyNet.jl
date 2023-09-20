import Base: ==

# Contains info for which species/taxa/etc. this lineage is an ancestor for
# while preserving coalescent order; e.g. [[3, 2], [1]]
# means that this is an ancestor of A, B, and C, and that BC must be
# together, coalescing before A
mutable struct Lineage
    # Fields
    lineage::Vector
    
    # Constructors
    Lineage(i::Real) = new([i])
    Lineage(l1::Lineage, l2::Lineage) = new(condense([l1.lineage, l2.lineage]))
    function Lineage(ls::Vector{Lineage}) 
        condensed = condense(ls)
        if typeof(condensed) <: Vector
            new(condensed)
        end
        new([condensed])
    end
    Lineage() = new([])
end

# Helper methods
ntaxa(lineage::Lineage) = length(lineage.lineage)
lineage(l::Lineage) = l.lineage
Base.coalesce(l1::Lineage, l2::Lineage) = Lineage(l1, l2)

# Condenses the `lineage` field as such:
# [[[[[1]]], [[[[3]]]]]] --> [1, 3]
function condense(ls)
    if typeof(ls) <: Lineage 
        return condense(ls.lineage)
    elseif typeof(ls) <: Vector && length(ls) == 1
        if typeof(ls[1]) <: Lineage
            return condense(ls[1])
        end
        return ls[1]
    elseif typeof(ls) <: Real
        return ls
    end
    return [condense(l) for l in ls]
end

# Returns all possible combinations of `ls` *in which all of the lineages
# end up coalescing*
#
# Order matters, so we can't just go from [1, 2, 3] -> [123]
# we need to instead go [1, 2, 3] -> [[1, 2], 3], [1, [2, 3]], [[1, 3], 2]
function Base.coalesce(ls::AbstractVector{Lineage})
    if length(ls) == 1 return ls end
    if length(ls) == 2 return [Lineage(ls[1], ls[2])] end

    choose2sets = partitions(ls, 2)
    retlist = []
    for set in choose2sets
        subsetsi = coalesce(set[1])
        subsetsj = coalesce(set[2])
        # subsetsi = ifelse(length(set[1]) > 2, coalesce(set[1]), set[1])
        # subsetsj = ifelse(length(set[2]) > 2, coalesce(set[2]), set[2])

        for opti in subsetsi
            for optj in subsetsj
                push!(retlist, Lineage(opti, optj))
            end
        end
    end
    return retlist
end

function Base.:(==)(l1::Lineage, l2::Lineage)
    if length(lineage(l1)) != length(lineage(l2)) return false end
    if length(lineage(l1)) == 1 return lineage(l1)[1] == lineage(l2)[1] end

    l2idxs = collect(1:length(lineage(l2)))
    for i in eachindex(lineage(l1))
        foundequality = -1
        for j in eachindex(lineage(l2))
            if lineage(l1)[i] == lineage(l2)[j]
                foundequality = j
                break
            end
        end
        if foundequality == -1 return false end
        deleteat!(l2idxs, foundequality)
    end

    return true
end

l1 = Lineage(Lineage(Lineage(1), Lineage(2)), Lineage(3))
l2 = Lineage(Lineage(3), Lineage(Lineage(1), Lineage(2))) 
l3 = Lineage(Lineage(1), Lineage(Lineage(3), Lineage(2)))
l1 == l2 || error("Lineage equality test failed.")
l1 != l3 || error("Lineage equality test failed.")
l2 != l3 || error("Lineage equality test failed.")