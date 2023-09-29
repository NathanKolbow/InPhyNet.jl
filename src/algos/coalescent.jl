# Gets all of the possible coalescent combinations for
# the lineages in `l.lineages`
#
# `bl` is the length of the branch that the
# coalescence is happening within
#
# returns a tuple (A, B)
# A: list of all possible coalescent combos
# B: the probability associated with the corresponding entry of A
function getcoalescentcombos(l::LineageNode, bl::Real)
    parts = partitions(lineages(l))
    partq = Queue{Any}()
    for part in parts enqueue!(partq, part) end

    N = nlineages(l)

    lineagelist = Vector{LineageNode}()
    problist = Vector{Float64}()
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

        # Calculate the probability of this outcome `sets`
        prob = _calculatecoalescentprobability(N, length(sets), bl)
        push!(problist, prob)

        # If we didn't get caught in the above loop, that means this partition is good to be pushed
        templin = LineageNode()
        for lin in sets
            push!(lineages(templin), Lineage(lin)) 
        end
        push!(lineagelist, templin)
    end

    return (lineagelist, problist)
end