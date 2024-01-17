# Need to map out how the results will be stored first.
include("retic-metric.jl")

using CSV, DataFrames

function savePerfectResults(truenet::HybridNetwork, constraints::Vector{HybridNetwork}, esterrors::Vector{<:Real},
    gausserrors::Vector{<:Real}, constraintdiffs::Vector{<:Real}, nretics_est::Vector{<:Real})

    # Quick checks for bad input
    a = length(gausserrors)
    b = length(esterrors)
    c = length(constraintdiffs)
    d = length(nretics_est)
    a == b == c == d || error("Input vector lengths not identical: $a - $b - $c - $d")

    # Relevant variables
    nrows = length(gausserrors)
    output_path = getOutputFilepath(truenet)
    if !isfile(output_path) touch(output_path) end
    
    # Calculate some relevant data points
    nretics_inside, nretics_outside, nretics_duplicated = calculateReticData(truenet, constraints)
    constraint_sizes = string([c.numTaxa for c in constraints])
    constraint_sizes = constraint_sizes[2:(length(constraint_sizes)-1)]

    # Put together and save data
    df = DataFrame(
        numtaxa=repeat([truenet.numTaxa], nrows),
        nretics_true=repeat([truenet.numHybrids], nrows),
        nretics_est=nretics_est,
        nretics_inside=repeat([nretics_inside], nrows),
        nretics_outside=repeat([nretics_outside], nrows),
        nretics_duplicated=repeat([nretics_duplicated], nrows),
        constraint_sizes=constraint_sizes,
        estRFerror=esterrors,
        gauss_error=gausserrors,
        constraint_error_sum=constraintdiffs
    )
    CSV.write(output_path, df, append=true)
end