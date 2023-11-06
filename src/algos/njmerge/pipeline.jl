"""
## Steps:

1. Calculate distance matrix $D \in \mathbb{R}_{\ge0}^{n\text{x}n}$
2. Compute quick, rough species tree $T$
3. Find possible reticulation relationships $R$
4. Decompose taxa into $k$ disjoint subsets $\mathbf{S}=\{S_1,\dots,S_k\}$ informed by $T$ and $R$
5. Estimate $m$ networks $\mathbf{N}=\{N_1,\dots,N_m\}$ on each $S_i$
6. Merge the set $\mathbf{N}_\text{C}$ into supernetwork $N_\text{super}$ from constraints $\mathbf{N}_\text{C}$
7. If desired, use $N_\text{super}$ in place of $T$ and repeat

## Inputs:
- msa: multiple sequence alignment
- <letter>fxn: each is 1-to-1 associated with the above steps (Dfxn 1st step, Tfxn 2nd step, ...)
- repeat: associated w/ 7th step above; `repeat>0` means that we repeat this pipeline `repeat` times, else we don't repeat
"""
function runpipeline(msa, Dfxn::Fxn, Tfxn::Fxn, Rfxn::Fxn, Sfxn::Fxn, Nfxn::Fxn, mergefxn::Fxn, repeat::Integer)
    Dout = Tout = Rout = Sout = Nout = mergeout = nothing
    
    if repeat <= 0 repeat = 1 end

    D = Dfxn(msa=msa)
    roughtree = Tfxn(msa=msa, D=D)
    retics = Rfxn(msa=msa, D=D, T=roughtree)
    subsets = Sfxn(msa=msa, D=D, T=roughtree, R=retics)
    constraints = Nfxn(msa=msa, D=D, T=roughtree, R=retics, S=subsets)
    supernet = mergefxn(msa=msa, D=D, T=roughtree, R=retics, S=subsets, N=constraints)

    for i in 2:repeat
        subsets = Sfxn(msa=msa, D=D, T=supernet, R=retics)
        constraints = Nfxn(msa=msa, D=D, T=supernet, R=retics, S=subsets)
        supernet = mergefxn(msa=msa, D=D, T=supernet, R=retics, S=subsets, N=constraints)
    end
end