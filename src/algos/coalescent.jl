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
        prob = _calculatecoalescentprobability(N, length(part), bl)

        # We need to make sure we're counting things like [1, 2], [3, 4] twice
        # or [1, 2], [3, 4], [5, 6] six times b/c order of coalescence is not
        # taken into account in the above `prob` calculation
        totalcoalevents = sum([length(set) >= 2 for set in part])
        prob *= factorial(totalcoalevents)
        
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

"""
Tests the above probabilities for `N` input lineages
and branch length `bl`.
"""
function probatest(N::Int64, bl::Real)
    N > 0 || error("N>0 required.")
    bl = BigFloat(bl)
    lnode = LineageNode([Lineage(i) for i=1:N])
    _, probs = getcoalescentcombos(lnode, bl)

    abs(sum(probs) - BigFloat(1.)) < 1e-12 || error("sum(probs) = "*string(sum(probs)))
end