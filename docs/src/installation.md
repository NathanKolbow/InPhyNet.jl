# Installation

## Install Julia

Julia is a high-level interactive programming language with C-like runtime speeds. To download and install Julia, [follow the instructions here](https://julialang.org/downloads/). For a quick introduction to Julia, see [this article from Learn X in Y minutes](https://learnxinyminutes.com/julia/).

## Install InPhyNet

To install InPhyNet, open an interactive Julia session and type:

```julia
using Pkg
Pkg.add("InPhyNet")
```

The walkthrough also uses [SNaQ.jl](https://github.com/JuliaPhylo/SNaQ.jl) to infer the constraint networks and [PhyloPlots.jl](https://github.com/JuliaPhylo/PhyloPlots.jl) to plot the final network. To follow along, install them as well:

```julia
Pkg.add(["SNaQ", "PhyloPlots"])
```
