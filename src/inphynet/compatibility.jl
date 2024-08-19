using PhyloNetworks

"""
Heuristic for determining whether the set of networks `ns` are compatible.
"""
function are_compatible_heuristic(ns::AbstractVector{<:HybridNetwork})
    for i = 1:(length(ns) - 1)
        for j = 2:length(ns)
            if !are_compatible_heuristic(ns[i], ns[j]) return false end
        end
    end
    return true
end

"""
Heuristic for determining whether trees `t1` and `t2` are compatible.
In the case of networks, `t1` and `t2` are the major trees of some other networks.
"""
function are_compatible_heuristic(t1::HybridNetwork, t2::HybridNetwork)
    leaf_overlap = intersect([leaf.name for leaf in t1.leaf], [leaf.name for leaf in t2.leaf])
    
    if length(leaf_overlap) <= 3 return true end

    t1_prime = pruneTruthFromDecomp(t1, leaf_overlap)
    t2_prime = pruneTruthFromDecomp(t2, leaf_overlap)
    
    return hardwiredClusterDistance(t1_prime, t2_prime, false) == 0
end


function are_compatible_after_merge(ns::AbstractVector{HybridNetwork}, nodenamei::AbstractString, nodenamej::AbstractString)

    for i = 1:(length(ns) - 1)
        i_leaves = [leaf.name for leaf in ns[i].leaf]
        if length(i_leaves) <= 3 continue end
        i_matches = (nodenamei in i_leaves) + (nodenamej in i_leaves)
        if i_matches == 0 continue end

        # Merge the nodes
        ni_prime = deepcopy(ns[i])
        if i_matches == 1
            nodej_match = findfirst(n -> n == nodenamej, i_leaves)
            if nodej_match !== nothing
                ni_prime.leaf[nodej_match].name = nodenamei
                i_leaves[nodej_match] = nodenamei
            end
        end

        for j = (i+1):length(ns)
            j_leaves = [leaf.name for leaf in ns[j].leaf]
            if length(j_leaves) <= 3 continue end
            j_matches = (nodenamei in j_leaves) + (nodenamej in j_leaves)
            if j_matches == 0 continue end
            
            # "Merge" the nodes
            nj_prime = deepcopy(ns[j])
            if j_matches == 1
                nodej_match = findfirst(n -> n == nodenamej, j_leaves)
                if nodej_match !== nothing
                    nj_prime.leaf[nodej_match].name = nodenamei
                    j_leaves[nodej_match] = nodenamei
                end
            end

            leaf_overlap = intersect(i_leaves, j_leaves)
            if length(leaf_overlap) <= 3 continue end

            nodej_match = findfirst(n -> n == nodenamej, leaf_overlap)
            if nodej_match !== nothing deleteat!(leaf_overlap, nodej_match) end

            if hardwiredClusterDistance(
                pruneTruthFromDecomp(ni_prime, leaf_overlap),
                pruneTruthFromDecomp(nj_prime, leaf_overlap),
                false) > 0 return false end
        end
    end

    return true
end
are_compatible_after_merge(n1::HybridNetwork, n2::HybridNetwork, nodenamei::AbstractString, nodenamej::AbstractString) = are_compatible_after_merge([n1, n2], nodenamei, nodenamej)
