<!-- [![Stable-Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://nathankolbow.github.io/InPhyNet.jl/stable) -->
<!-- [![Dev-Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://nathankolbow.github.io/inphynet.jl/dev) -->

# Table of Contents <!-- omit in toc -->
- [Installation](#installation)
- [Using InPhyNet](#using-inphynet)
  - [Inferring constraint networks](#inferring-constraint-networks)
  - [Calculating the dissimilarity matrix](#calculating-the-dissimilarity-matrix)
  - [Merging networks](#merging-networks)

&nbsp;

# Installation

```julia
using Pkg
Pkg.add("https://github.com/NathanKolbow/InPhyNet.jl")
```

# Using InPhyNet

## Inferring constraint networks

Constraint networks can be inferred with any software, but one solution that can be done in Julia alongside InPhyNet is SNaQ, which can be found in the [PhyloNetworks Julia package](https://juliaphylo.github.io/PhyloNetworks.jl/stable/).

## Calculating the dissimilarity matrix

Any method can be used to compute a dissimilarity matrix, though average gene tree internode distance (AGID) was used in our original simulation experiments and worked very well. InPhyNet provides function for calculating this distance matrix directly.

```julia
using PhyloNetworks, InPhyNet

# 1. Load the gene trees
gts = readMultiTopology("est_gts.treefile")

# 2. Calculate the dissimilarity matrix
D, namelist = calculateAGID(gts)
```

## Merging networks

After inferring your constraint networks and calculating your pairwise dissimilarity matrix `D` along with a `namelist` coinciding with each index of `D`, you can run InPhyNet.

```julia
using InPhyNet

# 1. We calculate D as above, though any method can be used
gts = readMultiTopology("est_gts.treefile")
D, namelist = calculateAGID(gts)

# 2. Load constraint networks
constraints = readMultiTopology("est_nets.netfile")

# 3. Run InPhyNet
full_net = netnj(D, constraints, namelist)
```
