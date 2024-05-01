function run_sim {
    netid=$1
    rep=$2
    m=$3
    dmethod=AGIC
    nsim=1000
    all_use_outgroup=$4
    remove_after_reroot=$5

    echo "${netid} ${rep} ${m} ${dmethod} ${all_use_outgroup} ${remove_after_reroot} ${nsim}"
    julia --project=../.. -t8 ./perfect_subsets.jl ${netid} ${rep} ${m} ${dmethod} ${all_use_outgroup} ${remove_after_reroot} ${nsim}
}

# Restricted params at first so that we can get results to look at
Nparallel=16
i=0
for rep in $(seq 1 25); do
    for top in n50r2 n50r5 n100r5 n100r10 n200r10 n200r20; do
        for maxsubsetsize in 5 10 15 20 25 30; do
            ((i=i%Nparallel)); ((i++==0)) && wait
            run_sim ${top} ${rep} ${maxsubsetsize} false false &
        done
    done
done

Nparallel=16
i=0
for rep in $(seq 1 25); do
    for top in n500r25 n500r50 n1000r50 n1000r100; do
        for maxsubsetsize in 5 10 15 20 25 30; do
            ((i=i%Nparallel)); ((i++==0)) && wait
            run_sim ${top} ${rep} ${maxsubsetsize} false false &
        done
    done
done
