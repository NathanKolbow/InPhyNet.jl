# Runs perfect subset sims for all params pairings w/ no noise added
include("helpers/helpers.jl")
InPhyNet.TIEWARNING = true

function run_perfect_no_noise_sim(net_id::String, rep_num::Int64, m::Int64)

    try
        true_net, cs, D, namelist = loadPerfectData(net_id, rep_num, m, "AGIC")
        mnet = netnj(D, cs, namelist)
        savePerfectResult(
            true_net,
            cs,
            getNetDistances(true_net, mnet),
            get_error_without_missing_retics(true_net, mnet),
            hardwiredClusterDistance(majorTree(true_net), majorTree(mnet), false),
            0.,
            0.,
            mnet.numHybrids,
            rep_num,
            m
        )
    catch e
        if typeof(e) <: SolutionDNEError
            savePerfectResult(
                true_net,
                cs,
                -1.,
                -1.,
                -1.,
                0.,
                0.,
                -1.,
                rep_num,
                m
            )
        end
    end
end

mutable struct AtomicCounter{Int64}; @atomic iterspassed::Int64; end
for net_id in reverse(get_all_net_ids())
    for m in [5, 10, 15, 20, 25, 30]
        @info "$(net_id), m=$(m)"
        ac = AtomicCounter(0)

        Threads.@threads for rep_num in reverse(collect(1:100))
            # Run sim
            run_perfect_no_noise_sim(net_id, rep_num, m)

            # Print progress
            @atomic :sequentially_consistent ac.iterspassed += 1
            if Threads.threadid() == 1
                print("\r\t$(ac.iterspassed) / 100")
            end
        end
        println()
    end
end