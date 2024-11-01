

```julia
using InPhyNet, PhyloNetworks

constraints = readMultiTopology(...)
gts = readMultiTopology(...)
D, namelist = calculateAGID(gts)

full_net = inphynet(D, constraints, namelist)
```