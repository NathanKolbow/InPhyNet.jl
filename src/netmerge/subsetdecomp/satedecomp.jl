# Implements subset decomposition methods outlines in SATe-I and SATe-II
using PhyloNetworks, StatsBase


"""
Performs subset decomposition as outlines in SATe-I on `inittree`.
"""
function sateIdecomp(inittree::HybridNetwork, maxsize::Integer)
    workingset = Vector{HybridNetwork}([inittree])
    splittrees = Vector{HybridNetwork}()

    while length(workingset) != 0
        curr = workingset[1]
        deleteat!(workingset, 1)

        medge = getMidSplitEdge(curr)
        split1, split2 = splitAtEdge(curr, medge)
        
        if split1.numTaxa > maxsize push!(workingset, split1)
        else push!(splittrees, split1) end
        
        if split2.numTaxa > maxsize push!(workingset, split2)
        else push!(splittrees, split2) end
    end

    return splittrees
end


"""
Performs subset decomposition as outlined in SATe-II on `inittree`.
"""
function sateIIdecomp(inittree::HybridNetwork, maxsize::Integer)
    workingset = Vector{HybridNetwork}([inittree])
    splittrees = Vector{HybridNetwork}()

    while length(workingset) != 0
        curr = workingset[1]
        deleteat!(workingset, 1)

        lbranch = getLongestBranch(curr)
        split1, split2 = splitAtEdge(curr, lbranch)
        
        if split1.numTaxa > maxsize push!(workingset, split1)
        else push!(splittrees, split1) end
        
        if split2.numTaxa > maxsize push!(workingset, split2)
        else push!(splittrees, split2) end
    end

    return splittrees
end


"""
Gets the longest edge in `tre`.
"""
function getLongestBranch(tre::HybridNetwork)
    # Find longest branch
    lbranch = tre.edge[1]
    for edge in tre.edge
        if edge.length > lbranch.length lbranch = edge end
    end
    return lbranch
end


"""
Finds a branch in `tre` that, when pruned on, splits the tree roughly in half.
"""
function getMidSplitEdge(tre::HybridNetwork)
    minsplitedge = nothing
    minsplitdiff = Inf

    for edge in tre.edge
        splitdiff = abs(length(getLeavesUnderEdge(edge)) - (tre.numTaxa / 2))
        if splitdiff < minsplitdiff && edge != tre.node[tre.root]
            minsplitedge = edge
            minsplitdiff = splitdiff
        end
    end

    return minsplitedge
end


"""
Splits `tre` into two trees at its longest branch.
"""
function splitAtEdge(tre::HybridNetwork, edge::PhyloNetworks.Edge)
    childset = [n.name for n in getLeavesUnderEdge(edge)]

    # Prune
    split1 = pruneTruthFromDecomp(tre, [childset])[1]
    split2 = pruneTruthFromDecomp(tre, [setdiff([l.name for l in tre.leaf], childset)])[1]

    return split1, split2
end


"""
Gets all the leaf nodes underneath edge `edge`.
"""
function getLeavesUnderEdge(edge::PhyloNetworks.Edge)
    children = Vector{PhyloNetworks.Node}()
    workingset = [getchild(edge)]

    while length(workingset) != 0
        curr = workingset[1]
        deleteat!(workingset, 1)

        if curr.leaf push!(children, curr)
        else
            for child in getchildren(curr) push!(workingset, child) end
        end
    end

    return children
end