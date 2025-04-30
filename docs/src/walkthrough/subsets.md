# Choosing your subsets

For thorough analyses, careful subset selection is important. InPhyNet cannot infer instances of reticulate evolution between its input networks, so several taxa from clades where reticulate evolution is likely to have occurred should appear together in at least one subset. To this end, we recommend first performing analyses with software like [TINNiK](https://cran.r-project.org/web/packages/MSCquartets/vignettes/TINNIK.html) to identify these areas.

Here, we take a naive approach and simply perform [centroid edge decomposition](https://pmc.ncbi.nlm.nih.gov/articles/PMC6705769/#:~:text=The%20clustering%20used%20in,each%20subset.&text=from%20the%20previous%20iteration%29%2C,each%20subset.&text=the%20two%20parts%20have,each%20subset.&text=each%20subtree%20until%20there,each%20subset.) on a species tree inferred by InPhyNet without any constraint networks (this is equivalent to Neighbor-Joining). Centroid edge decomposition is a recursive process that splits a "guide tree" in half at a single edge such that the two halves are of equal size. Then, this process is done recursively on the resulting trees until each subtree contains at most some predefined number of taxa.

This requires us to specify a maximum subset size–for this walkthrough we will use $8$. In practice we recommend that you use the highest number such that you can still infer networks in a reasonable amount of time with your preferred network inference software. We use $8$ here so that all of the walkthrough can be completed quickly.

Additionally, we specify a minimum subset size of $5$ because that is the smallest set of data that SNaQ can infer a network on.

```julia
nj_tre = inphynet(D, namelist)
subsets = centroid_edge_decomposition(nj_tre, 5, 8)
```

# Using TINNiK


