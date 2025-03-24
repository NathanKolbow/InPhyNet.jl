using PhyloNetworks

"""
Heuristic for determining whether the set of networks `ns` are compatible.
"""
function are_compatible_after_merge(ns::AbstractVector{HybridNetwork}, nodenamei::AbstractString, nodenamej::AbstractString)


    # Check 1: if only 1 or fewer networks has both nodes, they are definitely compatible
    relevant_nets = Vector{HybridNetwork}()
    for net in ns
        idx_i = findfirst(l -> l.name == nodenamei, net.leaf)
        idx_j = findfirst(l -> l.name == nodenamej, net.leaf)
        has_i = idx_i !== nothing
        has_j = idx_j !== nothing

        if !has_i && !has_j continue end
        net_copy = deepcopy(net)

        if !has_i && has_j
            net_copy.leaf[idx_j].name = nodenamei
        elseif has_i && has_j
            PhyloNetworks.deleteleaf!(net_copy, net_copy.leaf[idx_j])
        end
        push!(relevant_nets, net_copy)
    end
    if length(relevant_nets) <= 1 return true end

    # Check 2: for each pair of post-merge networks with overlapping taxa,
    #          when pruned to only contain their overlapping taxa, are
    #          they topologically identical? If not, this merge does not work.
    for i = 1:(length(relevant_nets) - 1)
        neti = relevant_nets[i]

        for j = (i+1):length(relevant_nets)
            netj = relevant_nets[j]

            overlapping_taxa = intersect(tiplabels(neti), tiplabels(netj))
            if length(overlapping_taxa) <= 3 continue end
            neti_pruned = prune_network(neti, overlapping_taxa)
            netj_pruned = prune_network(netj, overlapping_taxa)

            if hardwiredclusterdistance(neti_pruned, netj_pruned, false) > 0 return false end
        end
    end

    return true
end
are_compatible_after_merge(n1::HybridNetwork, n2::HybridNetwork, nodenamei::AbstractString, nodenamej::AbstractString) = are_compatible_after_merge([n1, n2], nodenamei, nodenamej)
