######################################## README ########################################
# In this file, we experiment to see which of the following gives better results:
# - taxa can be considered siblings if they are siblings in any displayed tree
# - taxa are only siblings if they are siblings in the major tree
########################################################################################

include("../../helpers/helpers.jl")


for net_id in ["n50r2", "n100r10", "n200r10"]
    for replicatenum in 1:5
        for m = [15, 20, 25]
            for major_tree_only in [true, false]
                for force_unrooted in [true, false]
                    @info "Running 100 simulations for $(net_id) m = $(m), $(major_tree_only)"
                    truenet, constraints, D, namelist = loadPerfectData(net_id, replicatenum, m, "AGIC")
                    esterrors, gausserrors, constraintdiffs, nretics_est =
                        monophyleticRobustness(truenet, constraints, D, namelist, nsim=100, displayprogress=true, major_tree_only=major_tree_only)

                    # Remove constraint error runs
                    keep_idxs = esterrors .!= -2.
                    esterrors = esterrors[keep_idxs]
                    gausserrors = gausserrors[keep_idxs]
                    constraintdiffs = constraintdiffs[keep_idxs]
                    nretics_est = nretics_est[keep_idxs]

                    @info "Saving results"
                    CSV.write("major_tree_only.csv",
                        DataFrame(
                            net_id=repeat([net_id], length(esterrors)),
                            hwcd=esterrors,
                            sum_constraint_hwcd=constraintdiffs,
                            gauss_sd=gausserrors,
                            nretics_est=nretics_est,
                            major_tree_only=repeat([major_tree_only], length(esterrors)),
                            m=repeat([m], length(esterrors)),
                            force_unrooted=repeat([force_unrooted], length(esterrors))
                    ), append=true)
                end
            end
        end
    end
end