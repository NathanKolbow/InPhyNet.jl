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
    retics_to_remove = intersect(retics_to_remove, [h.name for h in true_hybs])
    
    for retic_name in retics_to_remove
        # Sometimes removing 1 hybrid will also result in another being removed,
        # so we need to make sure the findfirst result is valid
        retic_idx = findfirst([h.name for h in true_hybs] .== retic_name)
        if retic_idx === nothing continue end
    
        hybnode = true_hybs[retic_idx]
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


"""
`true_net` is a baseline network, `est_net` is a network estimated all the way up from
sequence data (typically). This fxn finds the subset of reticulations in `true_net` such
that `true_net` has only as many retics as `est_net` and the HWCD is minimized between
the two networks.
"""
function find_minimum_retic_subset_hwcd(true_net::HybridNetwork, est_net::HybridNetwork)
    # Make sure the nets are properly rooted
    try_outgroup_root(true_net)
    try_outgroup_root(est_net)

    # If they have the same number of retics, or `est_net` somehow has more retics
    # than `true_net`, return their HWCD
    if est_net.numHybrids >= true_net.numHybrids return hardwiredClusterDistance(true_net, est_net, true) end

    # Find the subset!
    true_hyb_names = [h.name for h in true_net.hybrid]
    hyb_combinations = combinations(true_hyb_names, true_net.numHybrids - est_net.numHybrids)
    hwcds = Array{Float64}(undef, length(hyb_combinations))
    Threads.@threads for (i, hyb_subset_names) in collect(enumerate(hyb_combinations))
        # 1. Copy the true net
        true_net_copy = readTopology(writeTopology(true_net))
        
        # 2. Find the set of hybrids in `hyb_subset_names`
        remove_hybnodes = []
        for to_remove_name in hyb_subset_names
            push!(remove_hybnodes, true_net_copy.hybrid[findfirst([h.name == to_remove_name for h in true_net_copy.hybrid])])
        end

        # 3. Remove the set of hybrids from `true_net_copy`
        for hyb in remove_hybnodes
            PhyloNetworks.deletehybridedge!(true_net_copy, getparentedgeminor(hyb))
        end

        # 4. Compare w/ hardwiredClusterDistance, save if new minimum
        try
            rootatnode!(true_net_copy, "OUTGROUP")
        catch
        end
        hwcds[i] = hardwiredClusterDistance(true_net_copy, est_net, true)
    end
    return minimum(hwcds)
end


function try_outgroup_root(net::HybridNetwork)
    try
        rootatnode!(net, "OUTGROUP")
    catch
    end
end






function collect_retry_data(netfile_name::String)
    data_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data/"
    netid, replicatenum, ngt, seq_len, ils_level, maxsubsetsize, dmethod = split(netfile_name, "_")[2:8]
    netid = String(netid)
    replicatenum = parse(Int64, replicatenum)
    ngt = parse(Int64, ngt)
    seq_len = parse(Int64, seq_len)
    maxsubsetsize = parse(Int64, maxsubsetsize)
    ils_level = String(ils_level)

    est_constraints = readMultiTopology(joinpath(data_dir, netfile_name))
    est_gts = readMultiTopology(joinpath(data_dir, "estgt_$(netid)_$(replicatenum)_$(ngt)_$(seq_len)_$(ils_level).treefile"))
    est_D, est_namelist = calculateAGIC(est_gts)

    return est_constraints, est_D, est_namelist
end