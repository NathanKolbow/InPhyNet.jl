using CSV, DataFrames, Graphs, Combinatorics, Plots

################################################################
# These functions are actually related to subset decomposition #
################################################################

function requiredsubsets(df, cutoff)
    minsubset = minimumuniquesubset(df, cutoff=cutoff)
    return kconnected(df, minsubset, k=4, cutoff=0.01)
end

function minimumuniquesubset(df; cutoff=0.01)
    df = deepcopy(df)
    N = sum(df[:,"HBp_T3"] .<= cutoff)
    subset = []
    namelist = names(df)[2:(findfirst(names(df) .== "12|34")-1)]

    while sum(df[:,"HBp_T3"] .<= cutoff) > 0
        reject = df[:,"HBp_T3"] .<= cutoff
        counts = [sum(df[:,name] .== 1 .&& reject) for name in namelist if name in names(df)]
        
        removeidx = argmax(counts)
        myname = namelist[removeidx]

        numremoved = sum(df[:,myname] .== 1 .&& reject)
        removerows = df[:,myname] .== 1
        df = df[.!removerows,:]

        push!(subset, namelist[removeidx])
        deleteat!(namelist, removeidx)
        deleteat!(counts, removeidx)
    end

    return subset
end

function kconnected(df, minsubset; k::Int64=4, cutoff=0.01)
    k > 1 && k <= 4 || error("k must be in [2, 4]")
    kconnected = []
    namelist = names(df)
    namelist = namelist[2:(findfirst(namelist .== "12|34")-1)]

    for name in minsubset
        tempdf = view(df, df[:,name] .== 1 .&& df[:,"HBp_T3"] .<= cutoff, :)
        counts = [sum(tempdf[:,partnername] .== 1) for partnername in namelist]
        topk = sortperm(counts, rev=true)[1:k]
        push!(kconnected, namelist[topk])
    end
    return kconnected
end

################################################################
#                                                              #
################################################################

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