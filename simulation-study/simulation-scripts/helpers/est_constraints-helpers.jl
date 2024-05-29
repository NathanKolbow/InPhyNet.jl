



"""
Returns `true` if the AIC stopping rule has been met, `false` otherwise.

AIC stopping rule: if the AIC did not improve on the last iteration, stop. Else, continue.
"""
function should_stop_AIC(neg_logliks::AbstractVector{<:Real})
    if length(neg_logliks) < 2
        return false
    end

    calc_aic(i::Int64, nll::Real) = 2 * (i - 1 + nll)
    aic_prev = calc_aic(length(neg_logliks) - 1, neg_logliks[length(neg_logliks) - 1])
    aic_new = calc_aic(length(neg_logliks), neg_logliks[length(neg_logliks)])

    # Stop if we stayed the same or got worse
    return aic_new >= aic_prev
end


function empirical_snaq(s_idx::Int64, nhybrids::Int64, start_tree::HybridNetwork, snaq_df::DataCF)
    filename_prefix = snaq_filename_prefix(s_idx, nhybrids)
    runtime_file = "$(filename_prefix).runtimes"

    # If already inferred, just read from file
    # If not, run snaq
    if isfile(runtime_file)
        log("\t\t\tconstraint already inferred", :green)
    else
        snaq_runtime = @elapsed silently() do
            snaq!(start_tree, snaq_df, filename=filename_prefix, hmax=nhybrids, seed=42)
        end
        open(runtime_file, "w+") do f
            write(f, "$(snaq_runtime)")
        end
    end

    _, neg_Pll, snaq_runtime = parse_snaq_output("$(filename_prefix).out")
    log("\t\t\t-Ploglik = $(round(neg_Pll, digits=2))")
    log_runtime(snaq_runtime)

    return neg_Pll
end


function parse_snaq_output(filepath::AbstractString)
    prefix = filepath[1:(length(filepath)-4)]
    subset_number = split(prefix, "_")[4]
    subset_number = parse(Int64, subset_number[7:length(subset_number)])
    nhybrids = parse(Int64, split(prefix, "_")[5][2])

    snaq_line = readlines(filepath)[1]
    best_net = readTopology("$(split(snaq_line, ";")[1]);")
    neg_Ploglik = parse(Float64, split(snaq_line, "=")[2])

    snaq_runtime = Inf
    if isfile("$(prefix).runtimes")
        snaq_runtime = parse(Float64, readlines("$(prefix).runtimes")[1])
    end
    return best_net, neg_Ploglik, snaq_runtime
end


function log_runtime(rt_seconds::Float64)
    minutes = Int64(round(rt_seconds / 60, digits = 0))
    hours = Int64(round(rt_seconds / 60 / 60, digits = 0))
    days = Int64(round(rt_seconds / 60 / 60 / 24, digits = 0))

    log("\t\t\ttook $(days)d $(hours)h $(minutes)m", :green)
end