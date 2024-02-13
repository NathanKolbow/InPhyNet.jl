# Run from simulation-study/simulation-scripts/
for top in n50r2 n50r5 n100r5 # n100r10 ...
do
    for maxsubsetsize in 5 10 15 20 25
    do
        for rep in $(seq 1 100)
        do
            julia --project=../.. -t16 ./perfect_subsets.jl ${top} ${rep} ${maxsubsetsize} "internode_count" 1000 &>> perfect_subsets.log
        done
    done
done