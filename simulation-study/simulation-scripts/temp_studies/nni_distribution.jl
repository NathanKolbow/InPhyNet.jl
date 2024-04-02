# Small investigation to see whether there is significant difference between
# e.g. [1, 1, 1, 1, 1, 1] vs. [0, 0, 0, 0, 0, 6] NNI move distributions
#
# gaussSd for ALL is 1.0
include("../helpers/helpers.jl")


nsim = 500
netid = "n100r5"
truenet, constraints, D, namelist = loadPerfectData(netid, 2, 20, "internode_count")

esterrors = zeros(nsim, 3, 5) .- 1.
majortreeRFs = zeros(nsim, 3, 5) .- 1.
constraintdiffs = zeros(nsim, 3, 5) .- 1.

Threads.@threads for sim_iter in 1:nsim
    doprint = (sim_iter % 10) == 0 || sim_iter == 1
    for (idx_a, nets_with_nni) in enumerate([1, 3, 5])
        for (idx_b, num_nni_moves) in enumerate([1, 5, 15, 20, 35])
            which_have_nni = sample(1:length(constraints), nets_with_nni, replace=false)
            nnimoves = sample(which_have_nni, num_nni_moves, replace=true)
            nnimoves = Vector{Int64}([sum(nnimoves .== i) for i=1:length(constraints)])

            esterrors[sim_iter, idx_a, idx_b], majortreeRFs[sim_iter, idx_a, idx_b], cdiffs, _, _ = 
                runRobustSim(truenet, constraints, D, namelist, 1., 1., nnimoves)
            constraintdiffs[sim_iter, idx_a, idx_b] = sum(cdiffs, dims=1)[1]
        end
    end
end

# Save results
esterrors_n100 = esterrors
majortreeRFs_n100 = majortreeRFs
constraintdiffs_n100 = constraintdiffs


# n500r25
netid = "n500r25"
truenet, constraints, D, namelist = loadPerfectData(netid, 2, 20, "internode_count")

esterrors = zeros(nsim, 5, 5) .- 1.
majortreeRFs = zeros(nsim, 5, 5) .- 1.
constraintdiffs = zeros(nsim, 5, 5) .- 1.

Threads.@threads for sim_iter in 1:nsim
    doprint = (sim_iter % 10) == 0 || sim_iter == 1
    for (idx_a, nets_with_nni) in enumerate([1, 8, 15, 21])
        for (idx_b, num_nni_moves) in enumerate([1, 5, 15, 20, 35])
            which_have_nni = sample(1:length(constraints), nets_with_nni, replace=false)
            nnimoves = sample(which_have_nni, num_nni_moves, replace=true)
            nnimoves = Vector{Int64}([sum(nnimoves .== i) for i=1:length(constraints)])

            esterrors[sim_iter, idx_a, idx_b], majortreeRFs[sim_iter, idx_a, idx_b], cdiffs, _, _ = 
                runRobustSim(truenet, constraints, D, namelist, 1., 1., nnimoves)
            constraintdiffs[sim_iter, idx_a, idx_b] = sum(cdiffs, dims=1)[1]
        end
    end

end

# Save results
esterrors_n500 = esterrors
majortreeRFs_n500 = majortreeRFs
constraintdiffs_n500 = constraintdiffs


######## Save results to file ########
df = DataFrame(
    num_nets_with_moves = Vector{Int64}([]),
    num_total_moves =  Vector{Int64}([]),
    est_errors =  Vector{Float64}([]),
    constraint_diffs = Vector{Float64}([]),
    net_id = Vector{String}([])
)
for sim_iter in 1:nsim
    for (idx_b, num_nni_moves) in enumerate([1, 5, 15, 20, 35])
        # n100r5
        for (idx_a, nets_with_nni) in enumerate([1, 3, 5])
            push!(df, [nets_with_nni, num_nni_moves, 
                        esterrors_n100[sim_iter, idx_a, idx_b], 
                        constraintdiffs_n100[sim_iter, idx_a, idx_b],
                        "n100r5"])
        end
        # n500r25
        for (idx_a, nets_with_nni) in enumerate([1, 8, 15, 21])
            push!(df, [nets_with_nni, num_nni_moves,
                        esterrors_n500[sim_iter, idx_a, idx_b],
                        constraintdiffs_n500[sim_iter, idx_a, idx_b],
                        "n500r25"])
        end
    end
end
CSV.write("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/temp_studies/nni_distribution_data_2.csv", df)


# FIGURES IN `nni_distribution_figs.R`