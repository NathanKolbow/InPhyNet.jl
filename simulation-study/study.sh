#!/bin/bash

# Run from simulation-study/simulation-scripts/
for rep in $(seq 1 10)   # going to 10 right now just to get some preliminary results, we don't need every replicate yet...
do
    for maxsubsetsize in 5 10 25 30
    do
        for top in n50r2 n50r5 n100r5 n100r10 n200r10 n200r20 n500r25 n500r50 n1000r50 n1000r100
        do
            echo "julia --project=../.. -t4 ./perfect_subsets.jl ${top} ${rep} ${maxsubsetsize} "internode_count" 1000 2>> perfect_subsets.log"
            julia --project=../.. -t4 ./perfect_subsets.jl ${top} ${rep} ${maxsubsetsize} "internode_count" 1000 2>> perfect_subsets.log
        done
    done
done

# Small sample of runs to test new code changes...
for rep in $(seq 1 10)   # going to 10 right now just to get some preliminary results, we don't need every replicate yet...
do
    for maxsubsetsize in 5 10 15 20 25 30
    do
        echo "julia --project=../.. -t8 ./perfect_subsets.jl ${top} ${rep} ${maxsubsetsize} "internode_count" 1000 2>> perfect_subsets.log"
        julia --project=../.. -t8 ./perfect_subsets.jl ${top} ${rep} ${maxsubsetsize} "internode_count" 1000 2>> perfect_subsets.log
    done
done

# Small subset for testing out estimated_gts.jl
set -e
top=n50r2
for rep in 1
do
    for nloci in 100 1000
    do
        for seq_len in 100 1000
        do
            for ils_level in low high
            do
                echo "JULIA_DEBUG=Main julia --project=../.. -t8 ./estimated_gts.jl ${top} ${rep} ${nloci} ${seq_len} ${ils_level}"
                JULIA_DEBUG=Main julia --project=../.. -t8 ./estimated_gts.jl ${top} ${rep} ${nloci} ${seq_len} ${ils_level}
            done
        done
    done
done