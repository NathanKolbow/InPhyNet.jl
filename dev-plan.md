# Development Plan

## Algorithm Overview:

Input: MSA for $n$ species

Output: supernetwork with all $n$ species

1. Calculate distance matrix $D \in \mathbb{R}_{\ge0}^{n\text{x}n}$
2. Compute quick, rough species tree $T$
3. Find possible reticulation relationships $R$
4. Decompose taxa into $k$ disjoint subsets $\mathbf{S}=\{S_1,\dots,S_k\}$ informed by $T$ and $R$
5. Estimate $m$ networks $\mathbf{N}=\{N_1,\dots,N_m\}$ on each $S_i$
6. Merge the set $\mathbf{N}_\text{C}$ into supernetwork $N_\text{super}$ from constraints $\mathbf{N}_\text{C}$
7. If desired, use $N_\text{super}$ in place of $T$ and repeat

## Step details

| Step | Input | Output | Main method | Other prospective methods | Where is main method? |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| 1 | MSA | $D$ | AGID | other from Erin's paper, true top distance | unknown |
| 2 | MSA, [$D$] | $T$ | `PhyloNetworks.nj!` | IQTree, RAxML, true major tree | `julia` |
| 3 | MSA | $R$ | MSCQuartets | HyDe | R |
| 4 | $T$, $R$ | $\mathbf{S}$ | `novel` | None | DNE |
| 5 | $\mathbf{S}$, [MSA] | $\mathbf{N}_\text{C}$ | `SNaQ` | `PhyNest` | `julia` |
| 6 | $\mathbf{N}_\text{C}$, $D$ | $N_\text{super}$ | `novel` | None | DNE |
| 7 | $N_\text{super}$, $R$ | $\mathbf{S}$ | ----- | ----- | ----- |

## High level pseudo-code

Step 1: [this paper](https://academic.oup.com/sysbio/article/60/5/661/1644054?login=false)

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

Step 6:
- NJMerge but respect and maintain reticulations
- Choice: if taxa $A$ and $B$ are siblings only via reticulation, should they be allowed to merge via proposal? Kevin and I think probably, but could try playing with this.