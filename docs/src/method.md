# How InPhyNet works

InPhyNet builds one large phylogenetic network out of several small ones. It does
not infer networks itself. Instead, you infer small **constraint networks** on
subsets of taxa with any method you like (for example SNaQ, PhyNEST, PhyloNet or
NANUQ), and InPhyNet merges them into a single network, using a distance matrix on
all taxa to decide how they fit together. See the [InPhyNet paper](@ref citation)
for the details.

## Divide and conquer

A full InPhyNet analysis has four steps (the [Walkthrough](walkthrough/introduction.md)
goes through each one):

1. **Distance matrix.** Estimate a distance matrix ``D`` on all taxa, for example the
   average gene-tree internode distance ([`calculateAGID`](@ref)).
2. **Subsets.** Split the taxa into disjoint subsets, small enough for your network
   inference method.
3. **Constraint networks.** Infer one network on each subset, with any method.
4. **Merge.** Combine ``D`` and the constraint networks with [`inphynet`](@ref).

Network inference (step 3) is by far the slowest step, and it only ever runs on
small subsets. The subsets are independent, so they can be analyzed in parallel.
This is what lets InPhyNet scale to thousands of taxa.

Because the inference method is interchangeable, InPhyNet works with any type of
input data that your chosen method accepts, and it benefits directly from faster
or more accurate inference methods in the future.

## The merge step

The merge is a constrained version of **neighbor joining** (NJ):

1. Every taxon starts as its own cluster.
2. At each step, InPhyNet joins the pair of clusters with the best NJ criterion
   computed from ``D``, then updates ``D`` for the new cluster. Unlike plain NJ,
   it only considers pairs whose joining does not contradict any constraint
   network. By default (`use_heuristic=true`), it also avoids joins that would
   lead to a conflict later on.
3. When all taxa are joined, the reticulations from the constraint networks are
   added back, with their inheritance probabilities (γ).

So the relationships **within** each subset, including its reticulations, come from
its constraint network, and the way the subsets are connected **to each other**
comes from ``D``. With no constraint networks, `inphynet(D, namelist)` is plain NJ;
the walkthrough uses this to build a guide tree for choosing subsets.

With more than two constraint networks, the merge can occasionally reach a point
where no join is allowed. InPhyNet then automatically falls back to merging the
constraint networks two at a time, which always succeeds. See
[Pairwise InPhyNet](advanced/pairwise-usage.md).

## What InPhyNet can and cannot tell you

InPhyNet returns a level-1 network on all taxa, containing every reticulation from
the constraint networks with its γ.

It does **not**:

- **Find reticulations between subsets.** A reticulation can only appear in the
  merged network if it is in one of the constraint networks. Taxa from clades where
  you suspect reticulation should be in the same subset; see
  [Choosing your subsets](walkthrough/subsets.md).
- **Estimate branch lengths.** The merged network has no branch lengths.
- **Root the network.** The root of the merged network is arbitrary. If you have an
  outgroup, root the network on it with `rootatnode!`, as in
  [the walkthrough](@ref walkthrough_full_net).
- **Correct errors in the constraint networks.** Mistakes in a constraint network
  carry over into the merged network. In the paper's simulations, InPhyNet's
  accuracy tracked the accuracy of the inference method used for the constraint
  networks.

## Requirements on the input

- **Constraint networks** must be level-1 and binary: every internal node has
  exactly 3 edges, and the root has 2 (no polytomies). InPhyNet checks this and
  throws a `ConstraintError` otherwise.
- **Subsets** should be disjoint: each taxon appears in at most one constraint
  network.
- **The distance matrix** must include every taxon that appears in a constraint
  network, with names matching exactly. Taxa that appear in ``D`` but in no
  constraint network are allowed (InPhyNet warns about them) and are placed using
  ``D`` alone.
- **Multiple alleles per species:** if your gene trees have alleles as tips but
  your constraint networks have species, map alleles to species when computing
  ``D`` (see [Mapping alleles to species](walkthrough/d-matrix.md#Mapping-alleles-to-species)).

## Choosing the subset size

The merge itself is fast. Almost all the run time is in inferring the constraint
networks, which depends on the size of the largest subset, not on the total number
of taxa. Use the largest subset size that your inference method can handle in a
reasonable time: larger subsets mean fewer reticulations are missed between subsets.
