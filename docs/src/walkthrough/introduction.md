
# Introduction

This walkthrough will take you through the steps of generating each piece of data needed to run InPhyNet to infer a 50 taxa network rapidly and accurately. Data used at each step of this walkthrough is available via functions that are built-in to InPhyNet.jl, so feel free to follow along or use your own data instead!

---

### Requirements:
- Julia installed
- InPhyNet and SNaQ Julia packages installed

---

### Steps:



|| Step | Input | Output |
|-|------|-------|--------|
| 1 | [Estimate a pairwise distance matrix](d-matrix.md) | Input data $I^*$ | Matrix $D$ |
| 2 | [Separate your taxa into subsets](subsets.md) | $D$ | Subsets $\mathbf{S}=\{S_i\}_{i=1}^k$ |
| 3 | [Infer "constraint" networks on your subsets of taxa](constraints.md) | $I,\mathbf{S}$ | Networks $\mathbf{N}=\{N_i\}_{i=1}^k$ |
| 4 | [Put it together with InPhyNet](full_net.md) | $D,\mathbf{N}$ | Species network $\mathcal{N}$ |

$^*I$ can be any form of input data with which a distance matrix and semi-directed networks can be computed and inferred. Here, we utilize estimated gene trees.

