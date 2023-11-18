# make sure to activate project before running
cd("C:\\Users\\Nathan\\repos\\network-merging\\simulations\\n40h4")
include("../../src/main.jl")
softwarepath = "../../software"
ngt = 1000

# 1. Simulate gene trees under the true network
# 2. Simulate DNA sequences from the simulated gene trees
# 3. Estimate gene trees from these DNA sequences
if !isfile("./estgt.treefile")
    truenet = readTopology("n40h4.net")
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

# 5. Parse quartets
hybsubsets, treesubset = decomposeFromQuartets(hbptabfile, cutoff=0.01)

# 6. Estimate constraints with SNaQ
constraints = Array{HybridNetwork}(undef, length(hybsubsets))
for (j, hybsub) in enumerate(hybsubsets)
    if !isfile("./data/net$(j).networks")
        println("\n\n\t\tEstimating constraint $(j)\n")
        hybsub = hybsubsets[1]
        temptrees = Array{HybridNetwork}(undef, ngt)
        for (i, gt) in enumerate(estgts)
            temptrees[i] = deepcopy(gt)
            delleaves = []
            for leaf in temptrees[i].leaf
                if !(leaf.name in hybsub)
                    push!(delleaves, leaf)
                end
            end
            for leaf in delleaves
                PhyloNetworks.deleteLeaf!(temptrees[i], leaf)
            end
            temptrees[i].names = copy(hybsub)
        end
        q, t = countquartetsintrees(temptrees, showprogressbar=false)
        df = readTableCF(writeTableCF(q, t))
        
        # 8 SNaQ runs b/c that's how many processors we're using right now
        constraints[j] = snaq!(temptrees[1], df, hmax=Int64(ceil(length(hybsub) / 3)), filename="./data/net$(j)", runs=8)
    else
        cons = readlines("./data/net$(j).networks")
        constraints[j] = readTopology(split(cons[1], ", with -loglik")[1])
    end
end

# 7. Calculate distance matrix
D, namelist = calculateAGID(estgts)

# 8. Merge
mergednet = netnj!(D, constraints; namelist=namelist)