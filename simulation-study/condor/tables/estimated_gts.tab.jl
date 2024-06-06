using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/")

include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")
output_table = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/condor/tables/estimated_gts.tab"


open(output_table, "w+") do f
    lines_written = 0
    for net_id in ["n200r10", "n500r25"]
        for rep in 1:10
            for ngt in [100, 5000]
                for seq_len in [500, 1000]
                    for ils_level in ["low", "high"]
                        for m in [25]
                            if !estimated_sims_already_performed(net_id, rep, ngt, seq_len, ils_level, m)
                                write(f, "$(net_id),$(rep),$(ngt),$(seq_len),$(ils_level),$(m)\n")
                                lines_written += 1
                            end
                        end
                    end
                end
            end
        end
    end

    for net_id in ["n200r10"]
        for rep in 1:10
            for ngt in [100, 1000, 5000]
                for seq_len in [500, 1000]
                    for ils_level in ["low", "med", "high"]
                        for m in [15]
                            if !estimated_sims_already_performed(net_id, rep, ngt, seq_len, ils_level, m)
                                write(f, "$(net_id),$(rep),$(ngt),$(seq_len),$(ils_level),$(m)\n")
                                lines_written += 1
                            end
                        end
                    end
                end
            end
        end
    end
    @info "Wrote $(lines_written) to $(output_table)"
end