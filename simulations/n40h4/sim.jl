# make sure to activate project before running
cd("C:\\Users\\Nathan\\repos\\network-merging\\simulations\\n40h4")
include("../../src/main.jl")
softwarepath = "../../software"
ngt = 1000
truenet = readTopology("n40h4.net")

# 1. Simulate gene trees under the true network
# 2. Simulate DNA sequences from the simulated gene trees
# 3. Estimate gene trees from these DNA sequences
if !isfile("./estgt.treefile")
    simgts = simulatecoalescent(truenet, ngt, 1)

    Threads.@threads for i=1:length(simgts)
        gt = simgts[i]

        simgtfile = "./data/simgt_$(i).treefile"
        seqgenfile = "./data/seqgen_$(i).seq"
        writeTopology(gt, simgtfile)
        run(pipeline(`$(softwarepath)/seq-gen-windows.exe -s0.001 -n1 -f0.3,0.2,0.2,0.3 -mHKY -op -l1000 $simgtfile`, stdout=seqgenfile, stderr=devnull))

        try
            run(pipeline(`$(softwarepath)/iqtree-1.6.12-windows.exe -quiet -s $seqgenfile`, stdout=devnull, stderr=devnull))
        catch e
            run(pipeline(`$(softwarepath)/iqtree-1.6.12-windows.exe -quiet -s -redo $seqgenfile`, stdout=devnull, stderr=devnull))
        end
    end
end

# 4. Run MSCquartets on estimated gene trees
estgts = Array{HybridNetwork}(undef, ngt)
if isfile("./estgt.treefile")
    estgts = readMultiTopology("./estgt.treefile")
else
    for i=1:ngt estgts[i] = readTopology("./data/seqgen_$(i).seq.treefile") end
    writeMultiTopology(estgts, "./estgt.treefile")
end

estgtfile = "./estgt.treefile"
hbptabfile = "./estgt.hbptab"
if !isfile(hbptabfile)
    run(pipeline(`Rscript ../../src/algos/njmerge/mscquartets.R $(estgtfile) $(hbptabfile)`))
end

# 5. Calculate distance matrix
D, namelist = calculateAGID(estgts)

# 6. Parse quartets
hybsubsets, treesubset = decomposeFromQuartets(hbptabfile, cutoff=0.01)

# 7. Count quartets in ALL estimated gene trees
q, t = countquartetsintrees(estgts, showprogressbar=false)

# 8. Estimate constraints with SNaQ
constraints = Array{HybridNetwork}(undef, length(hybsubsets))
for (j, hybsub) in enumerate(hybsubsets)
    if !isfile("./data/net$(j).networks")
        # get only the quartets that are relevant to `hybsub`
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
    else
        cons = readlines("./data/net$(j).networks")
        constraints[j] = readTopology(split(cons[1], ", with -loglik")[1])
    end
end

# 9. Merge
mergednet = netnj!(D, constraints; namelist=namelist)
hardwiredClusterDistance(mergednet, truenet, false)


for n in mergednet.node
    if occursin("__internal", n.name)
        n.name = ""
    end
end
net = constraints[1]
for e in net.edge
    e.length = -1.
    if e.hybrid
        e.gamma = -1.
    end
end


# 10. Dig into results

## How accurate are the constraints?
## - RF dist, num hybrids
truec1 = 

## The # of hybrids is overestimated; what's our RF distance if we reduce to correct numbers?

using PhyloPlots
plot(truenet)
plot(mergednet)