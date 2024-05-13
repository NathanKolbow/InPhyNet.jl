# cd /mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/driver_scripts
# julia --project=../../.. -p10 -t10 ./perfect_subsets_driver.jl
# File for driving simulations, cleaner and easier than using many bash commands.

using Distributed

@everywhere include("../helpers/helpers.jl")
@everywhere InPhyNet.TIEWARNING = true

# pmap driver fxn
@everywhere function pmap_func(net_id::String, rep_num::Int64, m::Int64, dmethod::String, nsim::Int64=1000)
    n_already_performed = n_perfect_sims_already_performed(net_id, rep_num, m, false, false)
    nsim += 1 - n_already_performed
    if nsim <= 0
        return 0
    end

    # 1. gather ground truth network, constraint, distance matrix, and namelist
    truenet, constraints, D, namelist = loadPerfectData(net_id, rep_num, m, dmethod)

    seed = parse(Int64, "$(truenet.numTaxa)42$(truenet.numHybrids)42$(rep_num)")
    Random.seed!(seed)

    # 2. run robustness testing
    @info "Running sims w/ $(net_id) #$(rep_num), m=$(m)"
    esterrors, esterrors_without_missing_retics, majortreeRFs, gausserrors, constraintdiffs, nretics_est =
        monophyleticRobustness(truenet, constraints, D, namelist, nsim=nsim, displayprogress=false, do_no_noise_sim=(n_already_performed == 0))
    constraintdiffs = sum(constraintdiffs, dims=1)[1,:]

    # 3. filter out the results that had constraint errors
    keep_idxs = esterrors .!= -2.

    # 4. save results
    savePerfectResults(
        truenet,
        constraints,
        esterrors[keep_idxs],
        esterrors_without_missing_retics[keep_idxs],
        majortreeRFs[keep_idxs],
        gausserrors[keep_idxs],
        constraintdiffs[keep_idxs],
        nretics_est[keep_idxs],
        rep_num,
        m,
        false,
        false
    )

    return nsim
end

# values to pass to pmap_func
uq_net_ids = ["n50r2", "n50r5", "n200r10", "n200r20", "n100r10", "n100r5", "n500r25", "n500r50"]
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
shuffle_idxs = sample(1:length(pmap_net_ids), length(pmap_net_ids), replace=false)
pmap_net_ids = pmap_net_ids[shuffle_idxs]
pmap_rep_num = pmap_rep_num[shuffle_idxs]
pmap_m = pmap_m[shuffle_idxs]

# return value doesn't actually matter
pmap((net_id, rep_num, m) -> pmap_func(net_id, rep_num, m, "AGIC"), pmap_net_ids, pmap_rep_num, pmap_m, on_error = ex -> -1)