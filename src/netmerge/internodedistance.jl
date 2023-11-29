"""
    majorinternodedistance(N::HybridNetwork)

Calculates internode distances between all pairs of taxa in the major displayed tree
of network `N`. If `N` is a tree, then internode distances are as expected for a tree.
"""
function majorinternodedistance(N::HybridNetwork)
    return internodedistance(majorTree(N))
end


"""
    internodedistance(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}}=nothing)

Calculates internode distances between all pairs of taxa in network `N`.
"""
function internodedistance(N::HybridNetwork; namelist::Union{Nothing,<:AbstractVector{String}}=nothing)
    D = zeros(N.numTaxa, N.numTaxa)
    if namelist == nothing
        namelist = [l.name for l in N.leaf]
    end

    Ngraph = Graph(N)
    removeredundantedges!(Ngraph)
    nodelistidx = [findfirst([n.name == name for n in N.node]) for name in namelist]

    for i=1:(N.numTaxa-1)
        nodenumi = nodelistidx[i]
        nodei = N.node[nodenumi]
        for j=(i+1):N.numTaxa
            nodenumj = nodelistidx[j]
            nodej = N.node[nodenumj]

            D[i, j] = D[j, i] = length(a_star(Ngraph, nodenumi, nodenumj)) - 1
        end
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