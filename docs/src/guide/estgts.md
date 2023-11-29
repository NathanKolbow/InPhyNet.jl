Say we have estimated gene trees in the file `my-est-gts.treefile`. First, we use `MSCquartets` in `R` to check for potential hybrid relationships and to inform the subset decomposition algorithm later on.

> Note: `MSCquartets` is prohibitively slow for many, many taxa, so we are working on porting its methods to `julia` where they can run much faster.

```R
library(MSCquartets)
library(ape)

gtrees=read.tree(file="my-est-gts.treefile")
taxanames=taxonNames(gtrees)
QT=quartetTable(gtrees,taxanames)
RQT=quartetTableResolved(QT)
pTable=quartetTreeTestInd(RQT,"T3")
HolmBonferroni(pTable,"T3",cutoff)
write.csv(getHBpTable(treefile, cutoff=cutoff), "quartets.dat")
```

The results from `MSCquartets` are now in the file `quartets.dat`.

Now we can do everything else in `juila`.

```julia
using PhyloNetworks, NetMerge

mergednet = netnj("my-est-gts.treefile", "quartets.dat")
```

Alternatively, if you want to make manual adjustments to the default pipeline, you can break the pipeline down as below.

```julia
estgts = readMultiTopology("my-est-gts.treefile")

# Distance matrix
D, namelist = calculateAGID(estgts)

# Subset decomposition
hybsubsets, treesubset = decomposeFromQuartets("quartets.dat", cutoff=0.01)

# Estimating constraint networks with SNaQ
q, t = countquartetsintrees(estgts, showprogressbar=false)
startingtree = PhyloNetworks.nj!(deepcopy(D), deepcopy(namelist))

constraints = Array{HybridNetwork}(undef, length(hybsubsets))
for (j, hybsub) in enumerate(hybusbsets)
    # Filter out all quartets that are not exclusively composed of taxa in `hybsub`
    temptaxonnumbers = [i for i in 1:length(t) if t[i] in hybsub]
    tempq = view(q, [i for i in 1:length(q) if all([number in temptaxonnumbers for number in q[i].taxonnumber])])
    tempdf = readTableCF(writeTableCF(tempq, t))

    # Estimate the network
    constraints[j] = snaq!(startingtree, tempdf, hmax=Int64(ceil(length(hybsub) / 3)), runs=10)
end

mergednet = netnj(D, constraints; namelist=namelist)
```