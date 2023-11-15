
function majorinternodedistance(N::HybridNetwork)
    return internodedistance(majorTree(N))
end

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

function calculateAGID(Ns::AbstractVector{HybridNetwork})
    D, namelist = internodedistance(Ns[1])
    for j=2:length(Ns)
        D .+= internodedistance(Ns[j], namelist=namelist)[1]
    end
    return D ./ length(Ns), namelist
end