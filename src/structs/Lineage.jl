# Contains info for which species/taxa/etc. this lineage is an ancestor for
# while preserving coalescent order; e.g. [[3, 2], [1]]
# means that this is an ancestor of A, B, and C, and that BC must be
# together, coalescing before A
struct Lineage
    # Fields
    # lineage::Vector{Union{Int64, Vector}}
    lineage::Vector
    
    # Constructors
    Lineage(i::Real) = new([i])
    Lineage(l1::Lineage, l2::Lineage) = new([l1.lineage, l2.lineage])
    Lineage(ls::AbstractVector{Lineage}) = new(ls)
    Lineage() = new([])
end

# Helper methods
ntaxa(lineage::Lineage) = length(lineage.lineage)
Base.coalesce(l1::Lineage, l2::Lineage) = Lineage(l1, l2)

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
