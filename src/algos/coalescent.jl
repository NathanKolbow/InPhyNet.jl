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