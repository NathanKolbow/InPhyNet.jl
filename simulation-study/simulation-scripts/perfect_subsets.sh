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
Nparallel=10
i=0
for rep in $(seq 1 5); do
    for maxsubsetsize in 10 20 25; do
        for top in n50r2 n50r5 n100r5 n100r10 n200r10 n200r20; do
            for all_use_outgroup in true false; do
                if $all_use_outgroup; then
                    for remove_after_reroot in true false; do
                        ((i=i%Nparallel)); ((i++==0)) && wait
                        run_sim ${top} ${rep} ${maxsubsetsize} ${all_use_outgroup} ${remove_after_reroot} &
                    done
                else
                    ((i=i%Nparallel)); ((i++==0)) && wait
                    run_sim ${top} ${rep} ${maxsubsetsize} ${all_use_outgroup} ${remove_after_reroot} &
                fi
            done
        done
    done
done