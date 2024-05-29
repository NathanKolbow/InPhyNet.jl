# Loop through all sets of (net ID, replicate number, max subset size) adding 100 sims to each combo at a time indefinitely

include("helpers/helpers.jl")
InPhyNet.TIEWARNING = true  # disables the warning message when there are ties


function run_n_sims(netid::String, maxsubsetsize::Int64, replicatenum::Int64, all_have_outgroup::Bool=false, outgroup_removed_after_reroot::Bool=false, nsim::Int64=100)
    # 1. gather ground truth network, constraint, distance matrix, and namelist
    truenet, constraints, D, namelist = loadPerfectData(netid, replicatenum, maxsubsetsize, "AGIC", all_have_outgroup=all_have_outgroup, outgroup_removed_after_reroot=outgroup_removed_after_reroot)

    seed = parse(Int64, "$(truenet.numTaxa)42$(truenet.numHybrids)42$(replicatenum)")
    Random.seed!(seed)

    # 2. run robustness testing
    esterrors, esterrors_without_missing_retics, majortreeRFs, gausserrors, constraintdiffs, nretics_est =
        monophyleticRobustness(truenet, constraints, D, namelist, nsim=nsim, displayprogress=false, do_no_noise_sim=false)
    constraintdiffs = sum(constraintdiffs, dims=1)[1,:]

    # 3. filter out the results that had constraint errors
    keep_idxs = esterrors .!= -2.

    # 4. save results (thread safe)
    savePerfectResults(
        truenet,
        constraints,
        esterrors[keep_idxs],
        esterrors_without_missing_retics[keep_idxs],
        majortreeRFs[keep_idxs],
        gausserrors[keep_idxs],
        constraintdiffs[keep_idxs],
        nretics_est[keep_idxs],
        replicatenum,
        maxsubsetsize,
        all_have_outgroup,
        outgroup_removed_after_reroot
    )
    return length(keep_idxs)
end


# Loop forever
total_logged = AtomicCounter(0)

while true
    for rep in 1:100
        for m in [15, 20, 25, 5, 10, 30]
            for (all_have_outgroup, remove_outgroup) in ((false, false), (true, false), (true, true))
                # only a 1/10 chance to do (true, false) or (true, true) b/c this is going to just be 
                # left running and we are *far* less interested in those at this point
                if all_have_outgroup && rand() > 0.05 continue end

                Threads.@threads for net_id in get_all_net_ids()
                    try
                        n_logged = silently() do
                            run_n_sims(net_id, m, rep, all_have_outgroup, remove_outgroup, 100)
                        end
                        @atomic :sequentially_consistent total_logged.iterspassed += n_logged
                    catch e
                        println("")
                        @error "Error received with ($(net_id)#$(rep), m=$(m), ($(all_have_outgroup), $(remove_outgroup))"
                        println("")
                    end

                    if Threads.threadid() == 1
                        print("\rTotal logged: ")
                        printstyled("$(total_logged.iterspassed)", color=:cyan)
                    end
                end
            end
        end
    end
end