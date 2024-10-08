# TODO: change all return values to be NamedArrays instead of matrices and a vector of names

"""
    majorinternodedistance(N::HybridNetwork)

Calculates internode distances between all pairs of taxa in the major displayed tree
of network `N`. If `N` is a tree, then internode distances are as expected for a tree.
"""
function majorinternodedistance(N::HybridNetwork)
    return internodedistance(majorTree(N))
end


"""
    majorinternodecount(N::HybridNetwork)

Calculates major internode counts between all pairs of taxa in the major displayed tree
of network `N`. If `N` is a tree, then internode counts are as expected for a tree.
"""
function majorinternodecount(N::HybridNetwork)
    return internodecount(majorTree(N))
end


"""
    internodedistance(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}}=nothing)

Calculates internode distance between all pairs of taxa in network `N`.
"""
function internodedistance(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}}=nothing)
    N.numHybrids == 0 || throw(ErrorException("N must be tree-like."))

    D = pairwiseTaxonDistanceMatrix(N)
    if namelist === nothing
        idxs = sortperm(tipLabels(N))
        return D[idxs, idxs], tipLabels(N)[idxs]
    end

    # We want to re-order the idxs in `D` so that they match up with the provided namelist
    # So, if `namelist[j]` is "A", `map[j]` should return `i` where `tipLabels(N)[i]` is "A"
    tip_labels = tipLabels(N)
    label_map = Dict(taxon_name => j for (j, taxon_name) in enumerate(tip_labels))
    namelist_map = Dict(namelist_item => j for (j, namelist_item) in enumerate(namelist))

    idxfilter = Array{Int64}(undef, size(D, 1))
    for i = 1:size(D, 1)
        tip_name = tip_labels[i]
        idxfilter[label_map[tip_name]] = namelist_map[tip_name]
    end

    return view(D, idxfilter, idxfilter), namelist
    
end


"""
    internodecount(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}}=nothing)

Calculates internode counts between all pairs of taxa in network `N`.
WARNING: treats `N` as unrooted
"""
function internodecount(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}} = nothing)
    D = zeros(N.numTaxa, N.numTaxa)
    if namelist == nothing
        namelist = sort([l.name for l in N.leaf])
    end

    Ngraph = Graph(N, includeminoredges = false)
    removeredundantedges!(Ngraph, N, keeproot = false)
    # nodelistidx = [findfirst([n.name == name for n in N.node]) for name in namelist]
    nodelistidx = [findfirst(n -> n.name == name, N.node) for name in namelist]

    for i=1:(N.numTaxa-1)
        paths = bellman_ford_shortest_paths(Ngraph, nodelistidx[i])
        D[i, :] .= D[:, i] .= paths.dists[nodelistidx] .- 1
        D[i, i] = 0.
    end
    return D, namelist
end


"""
    calculate_pairwise_metric_across_networks(nets::AbstractVector{HybridNetwork}, net_metric::Function, summary_metric::Function)

Used as the base of `calculateAGIC`, `calculateAGID`, and `calculateMGIC`.

# Arguments
- `nets`: vector of `HybridNetwork` objects
- `net_metric`: metric calculated across each network - `net_metric(net::HybridNetwork, namelist::Vector{String}=namelist)`
                must be valid and return Float64
- `summary_metric`: metric that summarized values from `net_metric`. e.g. `mean` or `median`
"""
# function calculate_pairwise_metric_across_networks(nets::AbstractVector{HybridNetwork}, net_metric::Function, summary_metric::Function)
#     all_names = Set()
#     for net in Ns
#         union!(all_names, [leaf.name for leaf in net.leaf])
#     end
#     n = length(all_names)
#     namelist = sort(collect(all_names))

#     D = Matrix{Vector{Float64}}(undef, n, n)
#     D .= []

#     for net in nets
#         iter_names = sort([leaf.name for leaf in net.leaf])
#         idx_filter = findall(x -> x in iter_names, namelist)
#         iter_D = view(D, idx_filter, idx_filter)

#         pairwise_vals = net_metric(net, namelist=iter_names)
#         for i = 1:size(pairwise_vals, 1)
#             for j = 1:size(pairwise_vals, 2)
#                 push!(iter_D[i, j], pairwise_vals[i, j])
#             end
#         end
#     end
    
#     summary_matrix = Matrix{Float64}(undef, n, n)
#     for i = 1:size(D, 1)
#         for j = 1:size(D, 2)
#             summary_matrix[i, j] = summary_metric(D[i, j])
#         end
#     end

#     return summary_matrix, namelist
# end


"""
    calculateAGID(Ns::AbstractVector{HybridNetwork})

Calculates the **A**verage **G**ene tree **I**nternode **D**istance for all pairs of taxa
across all networks in `Ns`.

# Returns

1. AGID matrix
2. names of taxa corresponding to rows/columns in the AGID matrix
"""
function calculateAGID(Ns::AbstractVector{HybridNetwork}; allow_missing_pairs::Bool=false, default_missing_value=Inf)
    return calculate_average_network_metric(Ns, internodedistance, allow_missing_pairs=allow_missing_pairs, default_missing_value=default_missing_value)
end


"""
    calculateAGIC(Ns::AbstractVector{HybridNetwork})

Calculates the **A**verage **G**ene tree **I**nternode **C**ounts for all pairs of taxa
across all networks in `Ns`. If `allow_missing_pairs` is false, an error will be thrown if
there are pairs of taxa that never appear in `Ns` together. Otherwise, entries for taxa that
never appear together will be set to `default_missing_value`.

# Returns

1. AGIC matrix
2. names of taxa corresponding to rows/columns in the AGID matrix
"""
function calculateAGIC(Ns::AbstractVector{HybridNetwork}; allow_missing_pairs::Bool=false, default_missing_value=Inf)
    return calculate_average_network_metric(Ns, internodecount, allow_missing_pairs=allow_missing_pairs, default_missing_value=default_missing_value)
end


"""
Calculates `pairwise_metric` on each pair of networks in `Ns`. Called by `calculateAGIC` and `calculateAGID`.
`pairwise_metric` must be a function that takes a `HybridNetwork` as its first positional input, must accept
    the named argument `namelist`, and must return a 2-Tuple of (A, B), where A is the pairwise metric
    applied to all pairs of taxa in the supplied network, and B is the namelist used or generated by the fxn.
"""
function calculate_average_network_metric(Ns::AbstractVector{HybridNetwork}, pairwise_metric::Function; allow_missing_pairs::Bool=false, default_missing_value=Inf)
    all_names = Set()
    for net in Ns
        union!(all_names, tipLabels(net))
    end
    n = length(all_names)
    namelist = Vector{String}(sort(collect(all_names)))

    pair_appearance_count = zeros(UInt64, n, n)
    D = zeros(Float64, n, n)

    for j=1:length(Ns)
        iter_names = sort(tipLabels(Ns[j]))
        idx_filter = findall(i -> namelist[i] in iter_names, 1:length(namelist))
        iter_D = view(D, idx_filter, idx_filter)
        iter_count = view(pair_appearance_count, idx_filter, idx_filter)

        iter_D .+= pairwise_metric(Ns[j], namelist=iter_names)[1]
        iter_count .+= 1
    end

    _never_together = findall(e -> e == 0, pair_appearance_count)
    if length(_never_together) == 0
        return D ./ pair_appearance_count, namelist
    elseif !allow_missing_pairs
        throw(ArgumentError("$(Int(length(_never_together)/2)) pairs of taxa appear in 0 networks together."))
    else
        # Set to 1 to avoid the division error
        pair_appearance_count[_never_together] .= 1
        D = D ./ pair_appearance_count

        # Set default missing value
        D[_never_together] .= default_missing_value
        @warn "$(length(_never_together)) pairs of taxa appear in 0 networks together."
        return D, namelist
    end
end


function calculateMGIC(Ns::AbstractVector{HybridNetwork})
    D_one, namelist_one = internodecount(Ns[1])
    D_namelist_arr = Array{Matrix{Float64}}(undef, length(Ns))
    D_namelist_arr[1] = D_one
    n = size(D_one)[1]

    @info "A"
    for i=2:length(Ns)
        D_namelist_arr[i] = internodecount(Ns[i], namelist=namelist_one)[1]
    end

    @info "B"
    out_D = zeros(n, n)
    for i=1:(n-1)
        for j=(i+1):n
            vals = Array{Float64}(undef, length(Ns))
            for (k, D) in enumerate(D_namelist_arr)
                vals[k] = D[i, j]
            end
            out_D[i, j] = out_D[j, i] = median(vals)
        end
    end
    return out_D, namelist_one
end