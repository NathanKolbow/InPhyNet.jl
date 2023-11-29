
function pipelineSNaQMSCquartets(estgts::AbstractVector{HybridNetwork}, mscqfile::AbstractString; cutoff::Float64=0.01)
    error("Still need tree estimation method...")

    hybsubsets, treesubset = decomposeFromQuartets(mscqfile, cutoff=cutoff)
    q, t = countquartetsintrees(estgts, showprogressbar=false)

    # estimate constraints
    constraints = Array{HybridNetwork}(undef, length(hybsubsets)+length(treesubset))
    for (j, hybsub) in enumerate(hybsubsets)
        temptaxonnumbers = [i for i in 1:length(t) if t[i] in hybsub]
        tempq = view(q, [i for i in 1:length(q) if all([number in temptaxonnumbers for number in q[i].taxonnumber])])
        tempdf = readTableCF(writeTableCF(tempq, t))
        
        startingtree = nothing
        for gt in estgts
            startingtree = deepcopy(gt)
            delleaves = []
            for leaf in startingtree.leaf
                if !(leaf.name in hybsub)
                    push!(delleaves, leaf)
                end
            end
            for leaf in delleaves
                PhyloNetworks.deleteLeaf!(startingtree, leaf)
            end
            startingtree = readTopology(writeTopology(startingtree))
            if length(startingtree.names) == length(hybsub) && all([name in hybsub for name in startingtree.names])
                break
            else
                startingtree = nothing
            end
        end
        if startingtree === nothing error("pruning failed.") end

        constraints[j] = snaq!(startingtree, tempdf, hmax=Int64(ceil(length(hybsub) / 3)), filename="./data/net$(j)", runs=8)
    end

    D, namelist = calculateAGID(estgts)
    return netnj!(D, constraints; namelist=namelist)
end