# #!/bin/bash
# 
# Run this in a `tmux` session
# 
# cd /mnt/dv/wid/projects4/SolisLemus-network-merging/empirical-study/
# julia --project=/mnt/dv/wid/projects4/SolisLemus-network-merging/ -t6 ./est_constraints.jl

using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/")
cd("/mnt/dv/wid/projects4/SolisLemus-network-merging/")

using Dates
function log(msg::String, msg_color::Symbol=:black, timestamp_color::Symbol=:black)
    printstyled("[$(now())] ", color=timestamp_color)
    printstyled("$(msg)\n", color=msg_color)
end


# --------------------------------------------------------------------- #
m = 10
log("Using m=$(m)", :red)
# --------------------------------------------------------------------- #


log("Loading scripts", :cyan)
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/helpers/helpers.jl")
base_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/empirical-study/"
snaq_dir = joinpath(base_dir, "data", "snaq_data")
est_constraint_dir = joinpath(base_dir, "data", "est_constraint_data")
if !isdir(snaq_dir) mkdir(snaq_dir) end

# Load all the data
log("Loading gene trees", :cyan)
est_gt_file = joinpath(base_dir, "data", "Best.FAA.tre")
est_gts = readMultiTopology(est_gt_file)
log("Calculating D", :cyan)
est_D, est_namelist = calculateAGIC(est_gts)
log("Estimating NJ tree", :cyan)
nj_tre = silently() do
    nj(DataFrame(est_D, est_namelist))
end
log("Finding subsets", :cyan)
subsets = sateIdecomp(nj_tre, m)

# For each subset, infer networks with 0, 1, 2 reticulations
# Number of retics can be increased later for nets w/ continually increasing PSL
snaq_filename_prefix(subset_idx, nhyb) = joinpath(snaq_dir, "snaq_m$(m)_subset$(subset_idx)_h$(nhyb)")
log("Entering SNaQ loop - $(length(subsets) * 3) total networks to infer")
for (s_idx, subset) in enumerate(subsets)
    log("Subset $(s_idx):")
    log("\tlength(subset) = $(length(subset))")
    iter_gts = Vector{HybridNetwork}([])
    for (i, gt) in enumerate(est_gts)
        # Gene tree must have at least 4 of the taxa in `subset` to be included here
        taxa_in_gt = [leaf.name for leaf in gt.leaf]
        taxa_in_gt_and_subset = intersect(taxa_in_gt, subset)
        
        if length(taxa_in_gt_and_subset) >= 4
            push!(iter_gts, readTopology(writeTopology(pruneTruthFromDecomp(gt, subset))))
        end
    end

    # Get the starting tree
    log("\tcalculating initial tree", :cyan)
    iter_D, iter_namelist = calculateAGIC(iter_gts)
    start_tree = silently() do
        nj(DataFrame(iter_D, iter_namelist))
    end

    # Gather data for SNaQ
    log("\tgathering SNaQ data", :cyan)
    q, t = countquartetsintrees(iter_gts, showprogressbar=false)
    snaq_df = silently() do
        readTableCF(writeTableCF(q, t))
    end

    log("\tinferring constraints:", :cyan)
    for nhybrids=0:2
        log("\t\thmax = $(nhybrids) (network $((s_idx-1)*3+nhybrids+1)/$(length(subsets)*3))", :cyan)
        filename_prefix = snaq_filename_prefix(s_idx, nhybrids)
        runtime_file = "$(filename_prefix).runtimes"
        if isfile(runtime_file)
            log("\t\t\tconstraint already inferred", :green)
            continue
        end

        snaq_runtime = @elapsed silently() do
            snaq!(start_tree, snaq_df, filename=filename_prefix, hmax=nhybrids, seed=42)
        end
        open(runtime_file, "w+") do f
            write(f, "$(snaq_runtime)")
        end
        log("\t\t\ttook $(round(snaq_runtime / 60 / 60, digits=2)) hours ($(round(snaq_runtime / 60 / 60 / 24, digits=2)) days)", :green)
    end
end