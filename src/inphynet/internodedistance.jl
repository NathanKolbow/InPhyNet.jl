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
    D = zeros(N.numTaxa, N.numTaxa)
    if namelist == nothing
        namelist = Array{String}(undef, N.numTaxa)
        for i=1:N.numTaxa namelist[i] = N.leaf[i].name end
    end

    Ngraph, Nweights = Graph(N, withweights = true)
    nodelistidx = [findfirst([n.name == name for n in N.node]) for name in namelist]

    for i=1:(N.numTaxa-1)
        paths = bellman_ford_shortest_paths(Ngraph, nodelistidx[i], Nweights)
        D[i, :] .= D[:, i] .= paths.dists[nodelistidx]
        D[i, i] = 0.

        # for j=(i+1):(net.numTaxa)
        #     D[i, j] = D[j, i] = paths.dists[nodelistidx[j]] - 1
        # end
    end
    return D, namelist
end


"""
    internodecount(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}}=nothing)

Calculates internode counts between all pairs of taxa in network `N`.
"""
function internodecount(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}}=nothing)
    D = zeros(N.numTaxa, N.numTaxa)
    if namelist == nothing
        namelist = [l.name for l in N.leaf]
    end

    Ngraph = Graph(N)
    removeredundantedges!(Ngraph)
    nodelistidx = [findfirst([n.name == name for n in N.node]) for name in namelist]

    for i=1:(N.numTaxa-1)
        paths = bellman_ford_shortest_paths(Ngraph, nodelistidx[i])
        D[i, :] .= D[:, i] .= paths.dists[nodelistidx] .- 1
        D[i, i] = 0.
    end
    return D, namelist
end


"""
    calculateAGID(Ns::AbstractVector{HybridNetwork})

Calculates the **A**verage **G**ene tree **I**nternode **D**istance for all pairs of taxa
across all networks in `Ns`.

# Returns

1. AGID matrix
2. names of taxa corresponding to rows/columns in the AGID matrix
"""
function calculateAGID(Ns::AbstractVector{HybridNetwork})
    D, namelist = internodedistance(Ns[1])
    for j=2:length(Ns)
        D .+= internodedistance(Ns[j], namelist=namelist)[1]
    end
    return D ./ length(Ns), namelist
end


"""
    calculateAGIC(Ns::AbstractVector{HybridNetwork})

Calculates the **A**verage **G**ene tree **I**nternode **C**ounts for all pairs of taxa
across all networks in `Ns`.

# Returns

1. AGIC matrix
2. names of taxa corresponding to rows/columns in the AGID matrix
"""
function calculateAGIC(Ns::AbstractVector{HybridNetwork})
    D, namelist = internodecount(Ns[1])
    for j=2:length(Ns)
        D .+= internodecount(Ns[j], namelist=namelist)[1]
    end
    return D ./ length(Ns), namelist
end