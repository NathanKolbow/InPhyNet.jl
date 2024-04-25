######################################## README ########################################
# 
# In this file, we experiment to see which of the following gives better results:
# - InPhyNet run w/ normal subset decomposition from SATe-I
# - InPhyNet run w/ the outgroup placed in EVERY subset, after SATe-I chooses subsets
########################################################################################

include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")



function run_experiment(netid::String, replicate_num::Int64, max_subset_size::Int64; nsim::Int64=1000)
    # 1. Load data
    true_net_orig, constraints_default_orig, D_orig, namelist_orig = loadPerfectData(netid, replicate_num, max_subset_size, "AGIC")
    subsets = [[leaf.name for leaf in c.leaf] for c in constraints_default_orig]
    for s in subsets if !("OUTGROUP" in s) push!(s, "OUTGROUP") end end
    constraints_outgroup_orig = pruneTruthFromDecomp(true_net_orig, subsets)

    net_hwcd_default = zeros(nsim) .- 2.
    net_hwcd_outgroup = zeros(nsim) .- 2.
    constraint_diffs_default = zeros(nsim) .- 2.
    constraint_diffs_outgroup = zeros(nsim) .- 2.
    gauss_sigmas = zeros(nsim) .- 2.
    gauss_error_level = Array{String}(undef, nsim)

    Threads.@threads for sim_iter = 1:nsim
        # 2. Copy data
        true_net = deepcopy(true_net_orig)
        constraints_default = deepcopy(constraints_default_orig)
        constraints_outgroup = deepcopy(constraints_outgroup_orig)
        D = deepcopy(D_orig)
        namelist = deepcopy(namelist_orig)

        # 3. Perturb D
        std0 = upperTriangStd(D)
        sig = -1.
        rand_num = rand()
        if rand_num < 0.33
            gauss_error_level[sim_iter] = "low"
            sig = std0 / 4
        elseif rand_num < 0.66
            gauss_error_level[sim_iter] = "med"
            sig = std0 / 2
        else
            gauss_error_level[sim_iter] = "high"
            sig = std0
        end
        addnoise!(D, Normal(sig, sig))
        gauss_sigmas[sim_iter] = sig

        # 4. Perturb constraints
        totalnnimoves = Int64(round(rand(Uniform(0, 2 * length(constraints_default)))))
        nnimoves = sample(1:length(constraints_default), totalnnimoves, replace=true)
        nnimoves = Vector{Int64}([sum(nnimoves .== i) for i = 1:length(constraints_default)])
        iter_constraintdiffs_default = zeros(length(constraints_default))
        iter_constraintdiffs_outgroup = zeros(length(constraints_default))
        if maximum(nnimoves) > 0
            for (i, (c_default, c_outgroup, nmoves)) in enumerate(zip(constraints_default, constraints_outgroup, nnimoves))
                for _=1:nmoves doRandomNNI!(c_default) end
                for _=1:nmoves doRandomNNI!(c_outgroup) end
                iter_constraintdiffs_default[i] = hardwiredClusterDistance(constraints_default_orig[i], constraints_default[i], false)
                iter_constraintdiffs_outgroup[i] = hardwiredClusterDistance(constraints_outgroup_orig[i], constraints_outgroup[i], false)
            end
        end
        constraint_diffs_default[sim_iter] = sum(iter_constraintdiffs_default)
        constraint_diffs_outgroup[sim_iter] = sum(iter_constraintdiffs_outgroup)

        # 5. Re-root constraints w/ outgroups
        for c in constraints_outgroup
            rootatnode!(c, "OUTGROUP")
        end

        # 6. Run InPhyNet w/ default subsets
        try
            mnet_default = netnj(D, constraints_default, namelist, supressunsampledwarning=true)
            rootatnode!(mnet_default, "OUTGROUP")
            net_hwcd_default[sim_iter] = hardwiredClusterDistance(mnet_default, true_net_orig, true)
        catch
            net_hwcd_default[sim_iter] = -1.
        end

        # 7. Run InPhyNet w/ new subsets
        try
            mnet_new = netnj(D, constraints_outgroup, namelist, supressunsampledwarning=true)
            rootatnode!(mnet_new, "OUTGROUP")
            net_hwcd_outgroup[sim_iter] = hardwiredClusterDistance(mnet_new, true_net_orig, true)
        catch
            net_hwcd_outgroup[sim_iter] = -1.
        end
    end

    df = DataFrame(
        netid = repeat([netid], nsim),
        replicate_num = repeat([replicate_num], nsim),
        max_subset_size = repeat([max_subset_size], nsim),
        net_hwcd_default = net_hwcd_default,
        net_hwcd_outgroup = net_hwcd_outgroup,
        constraint_diffs_default = constraint_diffs_default,
        constraint_diffs_outgroup = constraint_diffs_outgroup,
        gauss_sigmas = gauss_sigmas,
        gauss_error_level = gauss_error_level
    )
    return df
end

output_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/temp_studies/all_constraints_have_outgroups.csv"
for rep_num in 1:10
    for max_subset_size in [25, 15]
        for net_id in ["n100r5", "n200r10", "n50r2"]
            @info "$(net_id): m=$(max_subset_size), rep=$(rep_num)"
            res_df = run_experiment(net_id, rep_num, max_subset_size, nsim = 100)
            CSV.write(output_file, res_df, append = isfile(output_file))
        end
    end
end