# Need to map out how the results will be stored first.
include("misc-metrics.jl")

using CSV, DataFrames

const AV = AbstractVector
function savePerfectResults(truenet::HybridNetwork, constraints::AV{HybridNetwork}, esterrors::AV{<:Real}, majortreeRFs::AV{<:Real},
    gausserrors::AV{<:Real}, constraintdiffs::AV{<:Real}, nretics_est::AV{<:Real}, replicate_num::Int64,
    max_subset_size::Real)

    # Quick checks for bad input
    a = length(gausserrors)
    b = length(esterrors)
    c = length(constraintdiffs)
    d = length(nretics_est)
    a == b == c == d || error("Input vector lengths not identical: $a - $b - $c - $d")

    # Relevant variables
    nrows = length(gausserrors)
    output_path = getOutputFilepath(truenet)
    if !isfile(output_path) copy_csv_template(output_path) end
    
    # Calculate some relevant data points
    nretics_inside, nretics_outside, nretics_duplicated = calculateReticData(truenet, constraints)
    constraint_sizes = string([c.numTaxa for c in constraints])
    constraint_sizes = replace(constraint_sizes[2:(length(constraint_sizes)-1)], " " => "")

    # Put together and save data
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
        majortreeRF=majortreeRFs,
        gauss_error=gausserrors,
        constraint_error_sum=constraintdiffs,
        replicate_num=repeat([replicate_num], nrows)
    )
    CSV.write(output_path, df, append=true)
end


function save_estimated_gts_results(netid::String, true_network::HybridNetwork, replicatenum::Int64, nloci::Int64,
    ils_level::String, max_subset_size::Int64, distance_method::String, seq_len::Int64,
    est_network::HybridNetwork, est_constraints::Vector{HybridNetwork}, est_gts::Vector{HybridNetwork},
    est_constraint_runtimes::Vector{Float64})

    # Basic info to make plotting easier
    num_taxa = split(netid, "r")
    num_retics_true = num_taxa[2]
    num_taxa = split(num_taxa[1], "n")[2]
    num_retics_est = est_network.numHybrids

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
        @debug "Finished HWCD calculation for mnet"

        # Calculate S
        @debug "Optimizing branch lengths of estimated network"
        q, t = countquartetsintrees(est_gts)
        df = readTableCF(writeTableCF(q, t))
        avg_bl = get_avg_bl(true_network)

        for e in est_network.edge       # optBL! needs non-negative starting branch lengths
            if e.hybrid && !e.isMajor
                e.length = 0.
            else
                e.length = avg_bl
            end
        end
        PhyloNetworks.optBL!(est_network, df)

        @debug "Starting S calculation"
        error("Need to optimize branch lengths on estimated network (because they are empty) before calcing mpl")
        S = calculate_net_logpseudolik(est_network, df) / calculate_net_logpseudolik(true_network, df)
        @debug "Finished S calculation"
    end

    @debug "Starting HWCD calculation for constraints"
    est_constraint_hwcd_sum = 0.
    for est_constraint in est_constraints
        taxa_names = [leaf.name for leaf in est_constraint.name]
        true_constraint = pruneTruthFromDecomp(true_network, taxa_names)
        est_constraint_hwcd_sum += hardwiredClusterDistance(true_constraint, est_constraint, false)
    end
    @debug "Finished HWCD calculation for constraints"

    # Create a DF then write the results
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

        # inferred newick
        est_newick = [est_newick],

        # inference errors
        est_net_hwcd = [est_net_hwcd],
        est_constraint_hwcd_sum = [est_constraint_hwcd_sum],
        S = [S],

        # time taken
        time_seconds_parallel = [maximum(est_constraint_runtimes)],
        time_seconds_series = [sum(est_constraint_runtimes)]
    )
end