#!/bin/bash

# Run from simulation-study/simulation-scripts/
for rep in $(seq 1 4)   # going to 4 right now just to get some preliminary results, we don't need every replicate yet...
do
    for maxsubsetsize in 5 10 15 20 25
    do
        for top in n50r2 n50r5 # n100r5 n100r10 ...
        do
            echo "julia --project=../.. -t6 ./perfect_subsets.jl ${top} ${rep} ${maxsubsetsize} "internode_count" 1000 2>> perfect_subsets.log"
            julia --project=../.. -t6 ./perfect_subsets.jl ${top} ${rep} ${maxsubsetsize} "internode_count" 10000 2>> perfect_subsets.log
        done
    done
done