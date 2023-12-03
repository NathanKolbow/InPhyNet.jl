function robustNNI(truenet::HybridNetwork, constraints::Vector{HybridNetwork},
    nmoves::Int64; nsim::Int64=100)

    dists = Array{Int64}(undef, nsim)
    constraintdiffs = Array{Int64}(undef, 4, nsim)
    

    # Threads.@threads for i=1:nsim
    #     truenet, constraints = loadTrueData("n40h4")
    #     for (j, c) in enumerate(constraints)
    #         origc = readTopology(writeTopology(c))
    #         for _=1:nmoves
    #             e = sample(c.edge, 1, replace=false)[1]
    #             nni!(c, e)
    #         end
    #         constraintdiffs[j, i] = hardwiredClusterDistance(origc, c, false)
    #     end

    #     mnet = runGroundTruthPipeline(truenet, constraints)
    #     dists[i] = hardwiredClusterDistance(truenet, mnet, false)
    # end

    Threads.@threads for i=1:nsim
        newconstraints = [readTopology(writeTopology(c)) for c in constraints]
        for (j, (c, newc)) in enumerate(zip(constraints, newconstraints))
            for _=1:nmoves
                e = sample(newc.edge, 1, replace=false)[1]
                nni!(newc, e)
            end
            constraintdiffs[j, i] = hardwiredClusterDistance(newc, c, false)
        end

        mnet = runGroundTruthPipeline(truenet, newconstraints)
        dists[i] = hardwiredClusterDistance(truenet, mnet, false)
    end
    return dists, constraintdiffs
end