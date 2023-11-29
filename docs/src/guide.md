```
using PhyloNetworks, NetMerge

estgts = readMultiTopology("my-est-gts.treefile")
mergednet = netmerge(estgts)
```