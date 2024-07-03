# Need to map out how the results will be stored first.
include("misc-metrics.jl")

using CSV, DataFrames, Pidfile

const AV = AbstractVector
function savePerfectResults(truenet::HybridNetwork, constraints::AV{HybridNetwork}, esterrors::AV{<:Real},
    esterrors_without_missing_retics::AV{<:Real}, majortreeRFs::AV{<:Real},
    gausserrors::AV{<:Real}, constraintdiffs::AV{<:Real}, nretics_est::AV{<:Real}, replicate_num::Int64,
    max_subset_size::Real, all_have_outgroups::Bool, outgroup_removed_after_reroot::Bool;
    level1::Bool=false)

    # Quick checks for bad input
    a = length(gausserrors)
    b = length(esterrors)
    c = length(constraintdiffs)
    d = length(nretics_est)
    a == b == c == d || error("Input vector lengths not identical: $a - $b - $c - $d")

    # Relevant variables
    nrows = length(gausserrors)
    output_path = level1 ? getOutputFilepathLevel1(truenet) : getOutputFilepath(truenet)
    if !isfile(output_path) copy_csv_template(output_path) end
    
    # Calculate some relevant data points
    nretics_inside, nretics_outside, nretics_duplicated = calculateReticData(truenet, constraints)
    constraint_sizes = string([c.numTaxa for c in constraints])
    constraint_sizes = replace(constraint_sizes[2:(length(constraint_sizes)-1)], " " => "")

    # Put together and save data
    lk = mkpidlock(output_path*".lk")
    try
        df = DataFrame(
            numtaxa=repeat([truenet.numTaxa - 1], nrows),
            nretics_true=repeat([truenet.numHybrids], nrows),
            nretics_est=nretics_est,
            nretics_inside=repeat([nretics_inside], nrows),
            nretics_outside=repeat([nretics_outside], nrows),
            nretics_duplicated=repeat([nretics_duplicated], nrows),
            constraint_sizes=constraint_sizes,
            max_subset_size=max_subset_size,
            estRFerror=esterrors,
            esterror_without_missing_retics=esterrors_without_missing_retics,
            majortreeRF=majortreeRFs,
            gauss_error=gausserrors,
            constraint_error_sum=constraintdiffs,
            replicate_num=repeat([replicate_num], nrows),
            all_have_outgroups=repeat([all_have_outgroups], nrows),
            outgroup_removed_after_reroot=repeat([outgroup_removed_after_reroot], nrows)
        )
        CSV.write(output_path, df, append=true)
    finally
        close(lk)
    end
end
savePerfectResult(true_net::HybridNetwork, constraints::AV{HybridNetwork},
    esterror::Real, esterror_without_missing_retics::Real, majortreeRF::Real,
    gausserror::Real, constraintdiff::Real, nretics_est::Real, rep_num::Int64,
    max_subset_size::Real, all_have_outgroup::Bool, outgroup_removed_after_reroot::Bool) = 
    savePerfectResults(true_net, constraints, [esterror],
        [esterror_without_missing_retics], [majortreeRF], [gausserror], [constraintdiff],
        [nretics_est], rep_num, max_subset_size, all_have_outgroup, outgroup_removed_after_reroot)


function save_estimated_gts_results(netid::String, true_network::HybridNetwork, replicatenum::Int64, nloci::Int64,
    ils_level::String, max_subset_size::Int64, distance_method::String, seq_len::Int64,
    est_network::HybridNetwork, est_constraints::Vector{HybridNetwork}, est_gts::Vector{HybridNetwork},
    est_constraint_runtimes::Vector{Float64}, inphynet_runtime::Float64)

    output_path = get_estimated_sim_output_filepath(true_network)
    if !isfile(output_path) copy_est_csv_template(output_path) end
    # Basic info to make plotting easier
    num_taxa = split(netid, "r")
    num_retics_true = num_taxa[2]
    num_taxa = split(num_taxa[1], "n")[2]
    num_retics_est = est_network.numHybrids

    # Retic data
    nretics_inside, nretics_outside, nretics_duplicated = calculateReticData(true_network, est_constraints)

    # Inference errors
    est_newick = ""
    est_net_hwcd = -1.
    S = -1.

    if est_network !== nothing
        @debug "Trying to re-root estimated network"
        try
            rootatnode!(est_network, "OUTGROUP")
        catch
        end
        @debug "Trying to re-root true network"
        try
            rootatnode!(true_network, "OUTGROUP")
        catch
        end

        @debug "Starting HWCD calculation for mnet"
        est_newick = writeTopology(est_network)
        est_net_hwcd = hardwiredClusterDistance(true_network, est_network, true)
        
        true_major = majorTree(true_network)
        est_major = majorTree(est_network)

        try
            rootatnode!(true_major, "OUTGROUP")
        catch
        end
        try
            rootatnode!(est_major, "OUTGROUP")
        catch
        end

        est_major_tree_hwcd = hardwiredClusterDistance(true_major, est_major, true)
    end

    @debug "Starting HWCD calculation for constraints"
    est_constraint_hwcd_sum = 0.
    for est_constraint in est_constraints
        taxa_names = [leaf.name for leaf in est_constraint.leaf]
        true_constraint = pruneTruthFromDecomp(true_network, taxa_names)
        est_constraint_hwcd_sum += hardwiredClusterDistance(true_constraint, est_constraint, false)
    end

    # @debug "Checking whether we should find the minimum retic subset HWCD"
    est_net_hwcd_only_identified_retics = "skipped"
    n_combinations = binomial(true_network.numHybrids, est_network.numHybrids)
    if n_combinations < 1e5
        est_net_hwcd_only_identified_retics = find_minimum_retic_subset_hwcd(true_network, est_network)
    else
        @debug "n_combinations = $(n_combinations) >= 1e5, skipping minimum retic subset HWCD step"
    end

    # Create a DF then write the results
    @debug "Creating DF"
    lk = mkpidlock(output_path*".lk")
    try
        df = DataFrame(
            # input param data
            net_id = [netid],
            replicate_num = [replicatenum],
            n_loci = [nloci],
            seq_len = [seq_len],
            ils_level = [ils_level],
            max_subset_size = [max_subset_size],
            distance_method = [distance_method],
    
            # basic topological info
            num_taxa = [num_taxa],
            num_retics_true = [num_retics_true],
            num_retics_est = [num_retics_est],
            num_retics_inside = [nretics_inside],
            num_retics_outside = [nretics_outside],
            num_retics_duplicated = [nretics_duplicated],
    
            # inferred newick
            est_newick = [est_newick],
    
            # inference errors
            est_net_hwcd = [est_net_hwcd],
            est_major_tree_hwcd = [est_major_tree_hwcd],
            est_net_hwcd_only_identified_retics = ["not implemented"],
            est_constraint_hwcd_sum = [est_constraint_hwcd_sum],
    
            # time taken
            time_seconds_parallel = [maximum(est_constraint_runtimes) + inphynet_runtime],
            time_seconds_series = [sum(est_constraint_runtimes) + inphynet_runtime]
        )
        @debug "Writing to file $(output_path)"
        CSV.write(output_path, df, append=true)
    finally
        close(lk)
    end
end


function estimated_sims_already_performed(netid::String, replicatenum::Int64, ngt::Int64, seq_len::Int64, ils_level::String, max_subset_size::Int64)

    output_path = get_estimated_sim_output_filepath(netid)
    if !isfile(output_path) return false end

    df = CSV.read(output_path, DataFrame)
    filt_df = filter(row -> row.net_id == netid && row.replicate_num == replicatenum &&
        row.n_loci == ngt && row.seq_len == seq_len && row.ils_level == ils_level && row.max_subset_size == max_subset_size, df)
    return nrow(filt_df) > 0
end


function n_perfect_sims_already_performed(netid::String, replicatenum::Int64, maxsubsetsize::Int64, all_have_outgroup::Bool, outgroup_removed_after_reroot::Bool)
    output_path = get_output_filepath(netid)
    if !isfile(output_path) return false end
    ntaxa = split(netid, "r")
    nretic = parse(Int64, ntaxa[2])
    ntaxa = parse(Int64, split(ntaxa[1], "n")[2])

    df = CSV.read(output_path, DataFrame)
    return nrow(filter(row -> row.replicate_num == replicatenum && row.max_subset_size == maxsubsetsize && row.all_have_outgroup == all_have_outgroup && row.outgroup_removed_after_reroot == outgroup_removed_after_reroot, df))
end


function no_noise_sim_already_performed(netid::String, replicatenum::Int64, maxsubsetsize::Int64, all_have_outgroup::Bool, outgroup_removed_after_reroot::Bool)
    output_path = get_output_filepath(netid)
    if !isfile(output_path) return false end
    ntaxa = split(netid, "r")
    nretic = parse(Int64, ntaxa[2])
    ntaxa = parse(Int64, split(ntaxa[1], "n")[2])

    df = CSV.read(output_path, DataFrame)
    return nrow(filter(row -> row.replicate_num == replicatenum && row.max_subset_size == maxsubsetsize && row.all_have_outgroup == all_have_outgroup && row.outgroup_removed_after_reroot == outgroup_removed_after_reroot && row.gauss_error == 0 && row.constraint_error_sum == 0, df)) > 0
end