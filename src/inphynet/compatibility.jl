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

    to_copy = Vector{Int64}([])
    ns_prime = Array{HybridNetwork}(undef, length(ns))
    n_matches = Array{Int64}(undef, length(ns))
    leaves = Array{Vector{String}}(undef, length(ns))
    
    # Check 1: if at most 1 network has nodes i and j, then that network is compatible
    #          with itself and we can just quit
    matching_ns = 0
    for (k, n) in enumerate(ns)
        leaves[k] = [leaf.name for leaf in n.leaf]
        n_matches[k] = (nodenamei in leaves[k]) + (nodenamej in leaves[k])
        
        if n_matches[k] > 0 && length(leaves[k]) > 3
            push!(to_copy, k)
            matching_ns += 1
        end
    end
    if matching_ns <= 1 return true end
    
    # Generate the pruned versions of the input networks
    for (k, n) in enumerate(ns)
        ns_prime[k] = deepcopy(n)
        
        if n_matches[k] == 1
            nodej_match = findfirst(n -> n == nodenamej, leaves[k])
            if nodej_match !== nothing
                ns_prime[k].leaf[nodej_match].name = nodenamei
                leaves[k][nodej_match] = nodenamei
            end
        elseif n_matches[k] == 2
            nodej_match = findfirst(n -> n == nodenamej, leaves[k])
            if nodej_match !== nothing
                deleteat!(leaves[k], nodej_match)
            end
        end
    end


    # Check 2: for each pair of post-merge networks with overlapping taxa,
    #          when pruned to only contain their overlapping taxa, are
    #          they topologically identical? If not, this merge does not work.
    #
    # TODO: we don't actually need to check every pair (i, j)...
    #       RF(t1, t2) = 0 implies that t1 and t2 are IDENTICAL,
    #       so we actually only need to check every pair (i, i+1)
    for i = 1:(length(ns) - 1)
        if !isassigned(ns_prime, i) continue end

        for j = (i+1):length(ns)
            if !isassigned(ns_prime, j) continue end

            leaf_overlap = intersect(leaves[i], leaves[j])
            if length(leaf_overlap) <= 3 continue end

            if hardwiredClusterDistance(
                pruneTruthFromDecomp(ns_prime[i], leaf_overlap),
                pruneTruthFromDecomp(ns_prime[j], leaf_overlap),
                false) > 0 return false
            end
        end
    end

    return true
end
are_compatible_after_merge(n1::HybridNetwork, n2::HybridNetwork, nodenamei::AbstractString, nodenamej::AbstractString) = are_compatible_after_merge([n1, n2], nodenamei, nodenamej)
