# MSCquartets implemented in Julia



# definitely have to read through things somewhat carefully, not just skim
function hbpgraph(hbpdf, cutoff)
    namelist = names(hbpdf)
    namelist = namelist[2:(findfirst(namelist .== "12|34")-1)]
    ntax = length(namelist)
    graph = SimpleGraph(ntax)
    hbpmat = Matrix(hbpdf)
    pvalcol = findfirst(names(hbpdf) .== "HBp_T3")

    for row in 1:size(hbpdf, 1)
        if hbpmat[row,pvalcol] <= cutoff
            connectedtaxa = findall(hbpmat[row,2:(ntax+1)] .== 1)
            for (taxai, taxaj) in combinations(connectedtaxa, 2)
                add_edge!(graph, taxai, taxaj)
            end
        end
    end
    return graph
end


# constraints[1] from mergealgo.jl; very connected, as expected
hbpdf = CSV.read("/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/10000-constraints1.hbptab", DataFrame)
graph = hbpgraph(hbpdf, 0.05)

connected_components(graph)
plot(graph)


#subnet1 = readTopology("(((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8)));")
#subnet2 = readTopology("(t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17)))))));")
net = readTopology("((((t1,(t2)#H1),((t3,#H1),t4)),((t5,(t6)#H2),((t7,#H2),t8))),(t9,(t10,((t11,#H3),(t12,(((t13,t14))#H3,(t15,(t16,t17))))))));")
for e in net.edge
    e.length = 1.
    if e.hybrid
        e.gamma = 0.5
    end
end
simfile = "/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/100-ex2.treefile"
writeMultiTopology(simulatecoalescent(net, 100, 1), simfile)
# then get hbptab from mscquartets.R
hbpdf = CSV.read("/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/100-ex2.hbptab", DataFrame)

plot(hbpgraph(hbpdf, 1/(100)))
plot(hbpgraph(hbpdf, 1/(100^2)))
plot(hbpgraph(hbpdf, 1/(100^3)))

connected_components(hbpgraph(hbpdf, 1/(100^2)))
connected_components(hbpgraph(hbpdf, 1/(100^3)))


hbpdf = CSV.read("/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/10000-ex2.hbptab", DataFrame)
plot(hbpgraph(hbpdf, 0.))

reqsubs = requiredsubsets(hbpdf, 0.01)