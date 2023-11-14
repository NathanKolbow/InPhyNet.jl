
function majorinternodedistance(N::HybridNetwork)
    D = zeros(N.numTaxa, N.numTaxa)
    names = [l.name for l in N.leaf]

    Ngraph = Graph(N, includeminoredges=false)
    removeredundantedges!(Ngraph)
    nodelistidx = [findfirst([n.name == name for n in N.node]) for name in names]

    for i=1:(N.numTaxa-1)
        nodenumi = nodelistidx[i]
        nodei = N.node[nodenumi]
        for j=(i+1):N.numTaxa
            nodenumj = nodelistidx[j]
            nodej = N.node[nodenumj]

            D[i, j] = D[j, i] = length(a_star(Ngraph, nodenumi, nodenumj)) - 1
        end
    end
    return D, names
end