# Gets all of the possible coalescent combinations for
# the lineages in `l.lineages`
#
# `bl` is the length of the branch that the
# coalescence is happening within
#
# returns a tuple (A, B)
# A: list of all possible coalescent combos
# B: the probability associated with the corresponding entry of A
"""
    getcoalescentcombos(l::LineageNode, bl::Real; complog::Union{CompDict,Nothing}=nothing)

Gets all possible coalescent combinations for the lineages in `l` and their associated 
probabilities, including ones where lineages in `l` do not coalesce.

# Arguments
- `l`: the relevant `LineageNode`
- `bl`: branch length that coalescing is happening in; used for probabilities
- `complog`: used to store computations so that the same work is not done multiple times
"""
function getcoalescentcombos(l::LineageNode, bl::Real, binoms::AbstractArray; complog::Union{CompDict,Nothing}=nothing)
    parts = partitions(lineages(l))
    N = nlineages(l)

    # We are going to be dividing by these values later
    counts = zeros(Int64, N)
    for part in parts counts[length(part)] += 1 end

    incompleteidxs = Vector{Int64}()
    lineagelist = Vector{LineageNode}()
    problist = Vector{BigFloat}()

    for (i, part) in enumerate(parts)
        # Calculate the probability of this outcome `part`
        push!(problist, _calculatetotalcoalescentprobability(N, length(part), bl, binoms, complog=complog) / counts[length(part)])

        templin = LineageNode()
        for lin in part
            push!(lineages(templin), Lineage(lin)) 
        end
        push!(lineagelist, templin)
    end

    return (lineagelist, problist)
end

function _findtotalcoaleventsinpart(part)
    totalevents = 1
    g2sets = 0
    for set in part
        totalevents *= _findtotalcoaleventsinpartrecur(set)
        g2sets += length(set) >= 2
    end
    println(max(totalevents, 1) * max(g2sets, 1))
    return max(totalevents, 1) * max(g2sets, 1)
end

function _findtotalcoaleventsinpartrecur(part)
    if typeof(part) <: Lineage return 1 end

    totalevents = 0
    nestedevents = 1
    for set in part
        if typeof(set) <: AbstractArray && length(set) >= 2
            totalevents += 1
            nestedevents *= max(_findtotalcoaleventsinpartrecur(set), 1)
        end
    end

    if totalevents == 0
        return binomial(length(part), 2)
    else
        return totalevents * nestedevents
    end
end


# bigl = [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
# _findtotalcoaleventsinpart(bigl)
# _findtotalcoaleventsinpart([[1, 2], [3]])
# _findtotalcoaleventsinpart([[1, 2], [3, 4]])
# _findtotalcoaleventsinpart([[1, 2], [3, 4], [5, 6]])
# _findtotalcoaleventsinpart([[1, 2], [3, 4], [5, 6, 7]])
# _findtotalcoaleventsinpart([[1, 2], [3, 4, 8], [5, 6, 7]])
# _findtotalcoaleventsinpart([[1, 2], 3, 4])
# _findtotalcoaleventsinpart([[1, [2, 3]], 3])
# _findtotalcoaleventsinpart([1, 2, [3, 4]])
# _findtotalcoaleventsinpart([1, 2, 3, 4])

# Next on the TODO list: run `probatoptest(4, 1.)` in arbitrary-proba.jl and look for incorrect things in the probs output

# Current state of errors: before, we realized the [1, 2], [3, 4] was getting counted once when it needed to get counted
# twice to account for the fact that either [1, 2] or [3, 4] could coalesce first. Now, we aren't double counting that in
# Base.coalesce in Lineage.jl by push!ing twice anymore, we're accounting for it in _findtotalcoaleventsinpart
# This isn't working as expected because some things still aren't getting counted the correct number of times
# but it is HORRIBLY UNCLEAR WHY THAT IS