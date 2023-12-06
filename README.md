[![Stable-Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://nathankolbow.github.io/network-merging/stable)
[![Dev-Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://nathankolbow.github.io/network-merging/dev)

# Development Plan <!-- omit in toc -->

## Table of contents <!-- omit in toc -->
- [Current Roadmap](#current-roadmap)
- [Algorithm overview](#algorithm-overview)
- [Step details](#step-details)
- [High level pseudo-code](#high-level-pseudo-code)
- [Repo Organization](#repo-organization)
- [Roadmap from 15 Nov 2023](#roadmap-from-15-nov-2023)
- [Roadmap from 06 Nov 2023](#roadmap-from-06-nov-2023)

&nbsp;

# Current Roadmap

Update 2 Dec 2023
- Robustness testing:
  - [ ] various number of NNI moves
  - [ ] adding noise to distance matrix
- [ ] Implement new subset decomposition methods discussed w/ Kevin and Claudia: estimate species tree w/ neighbor joining and do clade decomposition

# Algorithm overview

Input: MSA for $n$ species

Output: supernetwork with all $n$ species

1. Calculate distance matrix $D \in \mathbb{R}_{\ge0}^{n\text{x}n}$
2. Compute quick, rough species tree $T$
3. Find possible reticulation relationships $R$
4. Decompose taxa into $k$ disjoint subsets $\mathbf{S}=\{S_1,\dots,S_k\}$ informed by $T$ and $R$
5. Estimate $m$ networks $\mathbf{N}=\{N_1,\dots,N_m\}$ on each $S_i$
6. Merge the set $\mathbf{N}_\text{C}$ into supernetwork $N_\text{super}$ from constraints $\mathbf{N}_\text{C}$
7. If desired, use $N_\text{super}$ in place of $T$ and repeat

# Step details

| Step | Input | Output | Main method | Other prospective methods | Where is main method? |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| 1 | MSA | $D$ | AGID | other from Erin's paper, true top distance | unknown |
| 2 | MSA, [$D$] | $T$ | `PhyloNetworks.nj!` | IQTree, RAxML, true major tree | `julia` |
| 3 | MSA | $R$ | MSCQuartets | HyDe | R |
| 4 | $T$, $R$ | $\mathbf{S}$ | `novel` | None | DNE |
| 5 | $\mathbf{S}$, [MSA] | $\mathbf{N}_\text{C}$ | `SNaQ` | `PhyNest` | `julia` |
| 6 | $\mathbf{N}_\text{C}$, $D$ | $N_\text{super}$ | `novel` | None | DNE |
| 7 | $N_\text{super}$, $R$ | $\mathbf{S}$ | ----- | ----- | ----- |

# High level pseudo-code

Step 1:
- [How to calculate AGID](https://academic.oup.com/sysbio/article/60/5/661/1644054#115080675)

Step 4 (HyDe; triplets):
- $\mathbf{T}$ the set of all triplets from HyDe with $p$-value at most $p$; $T_{i,j}$ the $j$-th triplet where taxa $H_i$ is marked as a hybrid
- $H$ the set of all unique taxa that are hybrids in $\mathbf{T}$; $|H|=m$

```
function decompose(choices);
function score(x) = ?

T = triplets matrix; T[i,j] = T_{i,j}
choices = vector(m) = 1
best_choices = choices
best_score = score(choices)

while choices[m] <= length(T[m, .])
    increment(choices, 1)
    subsets = decompose(choices)
    if score(subsets) > best_score
        best_choices = choices
        best_score = score(subsets)

function increment(choices, idx)
    choices[idx] += 1
    if choices[idx] > length(T[idx, .])
        choices[idx] = 1
        if idx+1 <= m
            increment(choices, idx+1)
```

- `score` function properties:
  - all else the same, if a larger $p$-value is subbed in for a smaller one, the score increases
  - all else the same, a smaller $\max_{i}|S_i|$ gives a better score
- Three main options for the philosophy behind `score`:
    1. Take only the hybrid relationships w/ the best $p$-values and ignore the resultant decomposition subsets
    2. Look at all hybrid relationships w/ *acceptable* $p$-values and take the ones that lead to decomposition subsets that we expect to have the smallest overall runtime
    3. A hybrid of (1) and (2) that doesn't care exclusively about either $p$-values or runtimes, but instead takes a tradeoff somewhere in the middle

Step 4 (MSCQuartets; quartets):
- unsure at the moment, need to review how MSCQuartets works
- [MSCQuartets](https://academic.oup.com/bioinformatics/article/37/12/1766/5919582)
- [theory behind MSCQuartets](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7643865/)

Step 6:
- NJMerge but respect and maintain reticulations
- Choice: if taxa $A$ and $B$ are siblings only via reticulation, should they be allowed to merge via proposal? Kevin and I think probably, but could try playing with this.

```
D = distance matrix
H = dict[constraintedge] = (fromedge, toedge); dict that keeps track of where retics are linked
    - constraintedge is the edge object from the original constraint network
    - fromedge is the edge in mergednet where the reticulation will be going out from
    - toedge is the edge in mergednet where the reticulation will be going into

constraints = vector of constraint networks
mergednet = set of m disconnected nodes; merging them as we go

function initializeH(H, constraints)
    for net in constraints
        for retic in net
            H[retic] = (null, null)

function retainretics(i, j, ei, ej)
    if exists OUTGOING retic r directly above i
        H[r][1] = ei
    else if exists INCOMING retic r directly above i
        H[r][2] = ei

    if exists OUTGOING retic r directly above j
        H[r][1] = ej
    else if exists INCOMING retic r directly above j
        H[r][2] = ej

while not all nodes connected
    compute Q from D
    define x a placeholder for the new name of the subset containing merged taxa
    pairvalid = false
    for pairing (i, j) in sort(Q)
        for const_net in constraints
            if (i, j) in const_net AND (i, j) not siblings
                pairvalid = false
                skip (i, j)
            else
                pairvalid = true
        
    if !pairvalid error

    merge (i, j) in mergednet into x
        - ei the new edge coming up from i
        - ej the new edge coming up from j
    for const_net in constraints that contains i or j
        if const_net contains i AND j
            retainretics(i, j, ei, ej)
            merge (i, j) into x in const_net
        else if const_net contains i
            relabel i as x
        else
            relabel j as x

for retic in keys(H)
    draw edge from H[retic][1] to H[retic][2]
```

# Repo Organization

- `docs/` files for building the `julia` package documentation
- `simulations/`
  - `n<N>h<H>/` contains relevant simulation files for network with `<N>` taxa and `<H>` reticulations
  - `ground-truth/` simulation experiments where ground truth distances and constraints are used
  - `truth-with-NNI/` simulation experiments where ground truth distances and constraints are used, but constraint networks are permuted with random NNI moves
  - `all-estimated/` simulation experiments where everything is estimated from the ground up, starting with simulated DNA sequences
- `software/` software used for simulation experiments
- `src/` NetMerge package source code

# Roadmap from 15 Nov 2023

Notes: top 3 bullet points dropped as they are not a main focus.

- [ ] ~~Improve subset decomposition algorithm~~
- [ ] ~~Implement MSCquartets methods in pure `julia`~~
- [ ] ~~Test merging algorithm on larger, more complicated networks with true gene trees/constraints~~
- [X] Test merging algorithm in simulation with estimated gene trees & constraint networks

# Roadmap from 06 Nov 2023

- [ ] ~~Review MSCQuartets and fill in its pseudo-code implementation~~
- [X] Sketch high-level pseudo-code for our network version of NJMerge above
- [X] Sketch low-level pseudo-code for our network version of NJMerge above
- [X] Implement our network version of NJMerge in `julia`
- [X] Put together a high-level `pipeline` construct in `julia` where you can plug in the various options above
- [X] Implement AGID matrix calculation in `julia` (perhaps w/ gene trees as inputs)
- [X] Implement taxa decomposition schema in `julia`