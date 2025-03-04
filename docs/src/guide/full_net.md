

```julia
using InPhyNet, PhyloNetworks

constraints = readmultinewick(...)
gts = readmultinewick(...)
D, namelist = calculateAGID(gts)

full_net = inphynet(D, constraints, namelist)
```