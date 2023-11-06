# Network Merging Project Plan

## Goal

Design an algorithm that takes a set of $n$ estimated networks and merges them into a single network that includes all unique taxa from each of the $n$ networks.

## Background

Algorithms exist to estimate phylogenetic networks from sequences, estimated gene trees, quartet concordance factors, and more. These methods, though, are prohibitively slow and resource intensive for a large amount of taxa. We hope to bridge this computational gap by making it tangible to estimate many small networks and then merge them.

## Aims

1. **Create software to generate all parental trees belonging to a given network $N$ and their associated probabilities.**

With a set of input networks $\{N_i\}$, the individual network objects are difficult to work with. Trees, on the other hand, are easier to work with, and parental trees provide full information on a given network under the multi species network coalescent (MSNC) model.

2. **Develop an algorithm to reconcile differences in parental trees.**

With a network $N$ and its parental trees $P$, adding reticulations to some base tree $T$ until all of $P$ are obtainable from $T$ will result in a network that is close to or equal to $N$. Then, this algorithm can be given the set of parental trees that came from each $\{N_i\}$ to form a merged network. Similar methods have been attempted with displayed trees, but these methods have done poorly because they are not rooted in coalescent theory. Additionally, parental trees have more meaningful probabilities associated with them, so they offer a significant advantage over displayed trees.

3. **Incorporate branch length information into the estimated network.**

Reconciling differences in topologies provides a clear final network topology, but does not acknowledge differences in branch lengths among parental trees. So, a clear improvement upon Aim 2 would be to estimate branch lengths alongside topology.

## Significance

Evolutionary relationships between different species and taxa are important for many reasons, including cancer research, tracking changes in viruses and bacteria, and conservation efforts. Many evolutionary relationships are not adequately described in tree-like structures due to "horizontal events", such as migration, introgression, and hybridization. Researchers are often interested in the relationships between many different taxa, but network estimation methods are currently prohibitively slow when estimating a network with many taxa. A method that can merge smaller networks into one large network would bridge this computational gap, allowing for large-scale network inference.

## Innovation

No methods currently exist that take a set of smaller networks and use it to estimate a larger network.

## Deliverables

1. Open source `Julia` software for merging small networks into a large network.
2. Publications describing the methods used to merge networks in journals such as Bioinformatics
3. Conference presentations describing the methods herein

## Timeline

1. Software to generate parental trees and their associated probabilities should be written and tested within the Fall 2023 semester, before the end of November.
2. The algorithm to reconcile differences in parental trees should be written and its results should be evaluated before the end of the Spring 2024 semester.
3. Branch length information should be incorporated soon after Aim 2 is complete.

## Team Members

Aims 1 and 2 are coding projects with straightforward objective, so algorithm development could potentially be delegated to undergraduates. Aims 3 and 4 will definitely require brainstorming sessions with others, but it is unclear whether they will require other individuals spending significant amounts of time on the project.

## Potential Limitations and Alternatives

1. *Poor topology estimation accuracy*: If the developed method yields poor topological accuracy, then we need to develop an algorithm that more intelligently utilizes relationships between parental trees. E.g. quartets of parental trees could be weighted by their parental tree's probability and used with existing quartet based estimation methods.

## Future Work

1. This algorithm will likely not address reticulations high up in the final merged network, between taxa that came from different pre-estimated networks. A natural extension would be to figure out how to efficiently infer such reticulations.
2. Reconciling differences in parental trees is a good starting point, but it may not be enough on its own. Parental trees provide full information on their network under the MSNC, so a more involved algorithm could be developed that leverages these relationships.

## Other Miscellaneous Ideas

- Input networks $\{N_i\}$ with parental trees $\{P_{ij}\}$, get quartets from each $P_{ij}$ weighted by the probability of $P_{ij}$ under network $i$. Feed those quartets to `SNaQ`