# cd /mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/driver_scripts
# julia --project=../../.. -p10 -t10 ./perfect_subsets_driver.jl
# File for driving simulations, cleaner and easier than using many bash commands.

using Distributed

@everywhere include("../helpers/helpers.jl")
@everywhere InPhyNet.TIEWARNING = true

# pmap driver fxn
@everywhere function pmap_func(net_id::String, rep_num::Int64, m::Int64, dmethod::String, nsim::Int64=1000)
    if no_noise_sim_already_performed(net_id, rep_num, m, false, false)
        @info "\tSKIPPING $(net_id) ($(rep_num)), m=$m"
        return 0
    end
    @info "Running $(net_id) ($(rep_num)), m=$m"

    try
        # 1. gather ground truth network, constraint, distance matrix, and namelist
        truenet, constraints, D, namelist = loadPerfectData(net_id, rep_num, m, dmethod)
        mnet = netnj(D, constraints, namelist)
    
        # 2. Calculate error values
        try_outgroup_root(truenet)
        try_outgroup_root(mnet)
    
        esterror = getNetDistances(truenet, mnet)
        majortreeRF = hardwiredClusterDistance(majorTree(truenet), majorTree(mnet), false)
        gausserror = 0.
        constraintdiff = 0.
        esterror_without_missing_retics = get_error_without_missing_retics(truenet, mnet)
    
        nretics_est = mnet.numHybrids
    
        # 3. save results
        savePerfectResult(
            truenet,
            constraints,
            esterror,
            esterror_without_missing_retics,
            majortreeRF,
            gausserror,
            constraintdiff,
            nretics_est,
            rep_num,
            m,
            false,
            false
        )
    
        return nsim
    catch e
        @info "ERROR $(net_id) ($(rep_num)), m=$m"
        @info typeof(e)
        rethrow(e)
    end
end

# values to pass to pmap_func
uq_net_ids = ["n50r2", "n50r5", "n200r10", "n200r20", "n100r10", "n100r5", "n500r25", "n500r50", "n1000r100", "n1000r50"]
uq_rep_num = 1:100
uq_m = [5, 10, 15, 20, 25, 30]

# put these vals in long arrays
pmap_net_ids = []
pmap_rep_num = []
pmap_m = []
for n in uq_net_ids
    for r in uq_rep_num
        for m in uq_m
            push!(pmap_net_ids, n)
            push!(pmap_rep_num, r)
            push!(pmap_m, m)
        end
    end
end

# return value doesn't actually matter
pmap((net_id, rep_num, m) -> pmap_func(net_id, rep_num, m, "AGIC"), pmap_net_ids, pmap_rep_num, pmap_m, on_error = ex -> -1)