# Inferring constraint networks

Constraint networks can be inferred with any method of your choosing, such as [SNaQ](https://github.com/JuliaPhylo/SNaQ.jl), [PhyNEST](https://github.com/sungsik-kong/PhyNEST.jl), [PhyloNet](https://phylogenomics.rice.edu/html/phylonetTutorial.html), or [NANUQ (implemented in the MSCquartets R package)](https://www.rdocumentation.org/packages/MSCquartets/versions/2.0.1). We utilize SNaQ because it is also implemented in Julia, so the workflow can all be contained in a single, straight-forward script.

For each subset, we need to perform the following in order to infer a network with SNaQ:
1. Prune the estimated gene trees to only include taxa in the subset
2. Generate SNaQ's input data
3. Run SNaQ

When using SNaQ, it is best practice to (1) perform many runs for each inferred network (typically 10) and (2) infer networks with numbers of hybrids $h=0,1,2,3,...$ before performing model selection to determine the "best" network. For simplicity, here we only perform 2 runs for each inferred network and we only infer networks with $h=1$.

```julia
# First, make a folder to put the data in
isdir("snaq_data") || mkdir("snaq_data")

# We will store the inferred networks in this array
snaq_networks = Array{HybridNetwork}(undef, length(subsets))

for (i, subset) in enumerate(subsets)
    @info "Inferring network for subset $(i)/$(length(subsets))"

    # 1. Prune estimated gene trees
    subset_gts = prune_networks(est_gts, subset)

    # 2. Generate SNaQ's input data
    df = readtrees2CF(subset_gts, writeTab=false, writeSummary=false)

    # 3. Run SNaQ
    snaq_networks[i] = snaq!(subset_gts[1], df, hmax=1, filename="snaq_data/snaq_subset$(i)", runs=1, seed=42)
end
```

