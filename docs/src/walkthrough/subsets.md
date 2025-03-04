# Choosing your subsets

For thorough analyses, careful subset selection is important. InPhyNet cannot infer instances of reticulate evolution between its input networks, so several taxa from clades where reticulate evolution is likely to have occurred should appear together in at least one subset. To this end, we recommend first performing analyses with software like [TINNiK](https://cran.r-project.org/web/packages/MSCquartets/vignettes/TINNIK.html) to identify these areas.

Here, we take a naive approach and simply perform centroid edge decomposition on a species tree inferred by InPhyNet without any constraint networks (this is equivalent to Neighbor-Joining). This requires us to specify a maximum subset size–for this walkthrough we will use $8$. In practice we recommend that you use the highest number such that you can still infer networks in a reasonable amount of time with your preferred network inference software. We use $8$ here so that all of the walkthrough can be completed quickly.

Additionally, we specify a minimum subset size of $5$ because that is the smallest set of data that SNaQ can infer a network on.

```julia
nj_tre = inphynet(D, namelist)
subsets = centroid_edge_decomposition(nj_tre, 5, 8)
```

