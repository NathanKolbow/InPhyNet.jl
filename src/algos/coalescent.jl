# Gets all of the possible coalescent combinations for
# the lineages in `l.lineages`
#
# `bl` is the length of the branch that the
# coalescence is happening within
#
# returns a tuple (A, B)
# A: list of all possible coalescent combos
# B: the probability associated with the corresponding entry of A
function getcoalescentcombos(l::LineageNode, bl::Real; complog::Union{CompDict,Nothing}=nothing)
    parts = partitions(lineages(l))
    partq = Queue{Any}()
    for part in parts enqueue!(partq, part) end

    N = nlineages(l)

    lineagelist = Vector{LineageNode}()
    problist = Vector{BigFloat}()
    while !isempty(partq)
        part = dequeue!(partq)

        continuewhile = false
        for (i, set) in enumerate(part)
            if length(set) > 2
                # If this grouping has > 2 taxa then we need to break it up further
                outcomes = coalesce(set)
                for (j, outcome) in enumerate(outcomes)
                    duppart = ifelse(j == length(outcomes), part, deepcopy(part))
                    duppart[i] = [outcome]

                    # re-queue the new copies
                    enqueue!(partq, duppart)
                end

                # We've re-queued the updated objects, so move on in the queue
                continuewhile = true
                break
            end
        end

        if continuewhile continue end

        # Calculate the probability of this outcome `part`
        prob = _calculatecoalescentprobability(N, length(part), bl, complog=complog)

        # We need to make sure we're counting things like [1, 2], [3, 4] twice
        # or [1, 2], [3, 4], [5, 6] six times b/c order of coalescence is not
        # taken into account in the above `prob` calculation
        totalcoalevents = _findtotalcoaleventsinpart(part)
        prob *= totalcoalevents
        
        push!(problist, prob)

        # If we didn't get caught in the above loop, that means this partition is good to be pushed
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


bigl = [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
_findtotalcoaleventsinpart(bigl)
_findtotalcoaleventsinpart([[1, 2], [3]])
_findtotalcoaleventsinpart([[1, 2], [3, 4]])
_findtotalcoaleventsinpart([[1, 2], [3, 4], [5, 6]])
_findtotalcoaleventsinpart([[1, 2], [3, 4], [5, 6, 7]])
_findtotalcoaleventsinpart([[1, 2], [3, 4, 8], [5, 6, 7]])
_findtotalcoaleventsinpart([[1, 2], 3, 4])
_findtotalcoaleventsinpart([[1, [2, 3]], 3])
_findtotalcoaleventsinpart([1, 2, [3, 4]])
_findtotalcoaleventsinpart([1, 2, 3, 4])

# Next on the TODO list: run `probatoptest(4, 1.)` in arbitrary-proba.jl and look for incorrect things in the probs output

# Current state of errors: before, we realized the [1, 2], [3, 4] was getting counted once when it needed to get counted
# twice to account for the fact that either [1, 2] or [3, 4] could coalesce first. Now, we aren't double counting that in
# Base.coalesce in Lineage.jl by push!ing twice anymore, we're accounting for it in _findtotalcoaleventsinpart
# This isn't working as expected because some things still aren't getting counted the correct number of times
# but it is HORRIBLY UNCLEAR WHY THAT IS