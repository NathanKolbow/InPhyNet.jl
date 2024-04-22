"""
Calculates the `nretics_inside`, `nretics_outside`, and `nretics_duplicated` metrics.
"""
function calculateReticData(truenet::HybridNetwork, constraints::Vector{HybridNetwork})
    true_retic_names = [retic.name for retic in truenet.hybrid]
    constraint_retic_names = Set{String}()

    nretics_duplicated = 0
    for c in constraints
        # Log retics, counting duplicates
        c_reticnames = [retic.name for retic in c.hybrid]
        for retic in c_reticnames
            if retic in constraint_retic_names
                nretics_duplicated += 1
            else
                push!(constraint_retic_names, retic)
            end
        end
    end

    nretics_inside = length(constraint_retic_names)
    nretics_outside = truenet.numHybrids - nretics_inside

    return nretics_inside, nretics_outside, nretics_duplicated
end


"""
Calculates the network `net`'s pseudo-likelihood given gene trees `gts`.
Code taken from PhyloNetworks `pseudolik.jl`
"""
function calculate_net_logpseudolik(net::HybridNetwork, df::PhyloNetworks.DataCF)
    net0 = deepcopy(net)
    PhyloNetworks.parameters!(net0)
    Threads.@threads for q in df.quartet
        PhyloNetworks.extractQuartet!(net0, q)
        PhyloNetworks.calculateExpCFAll!(q.qnet)
    end

    return PhyloNetworks.logPseudoLik(df)
end


"""
Gets the HWCD error between `true_net` and `mnet` after removing reticulations in
`true_net` that do not appear in `mnet`. Hybrids must match in name for this
function to work properly.
"""
function get_error_without_missing_retics(true_net::HybridNetwork, mnet::HybridNetwork; try_root_at_node::String="OUTGROUP")
    true_net_copy = readTopology(writeTopology(true_net))
    true_hybs = true_net_copy.hybrid
    mnet_hybs = mnet.hybrid
    retics_to_remove = setdiff([h.name for h in true_hybs], [h.name for h in mnet_hybs])

    for retic_name in retics_to_remove
        hybnode = true_hybs[findfirst([h.name for h in true_hybs] .== retic_name)]
        PhyloNetworks.deletehybridedge!(true_net_copy, getparentedgeminor(hybnode))
    end

    try
        rootatnode!(true_net, "OUTGROUP")
    catch
    end
    try
        rootatnode!(mnet, "OUTGROUP")
    catch
    end
    
    return hardwiredClusterDistance(true_net_copy, mnet, true)
end