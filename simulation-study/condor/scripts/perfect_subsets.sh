#!/bin/bash
# $1: netid
# $2: replicate number
# $3: maximum subset size
# $4: distance method ("internode_count")
# $5: number of sims

netid=$1
replicate=$2
maxsubsetsize=$3
dmethod=$4
nsim=$5
nthreads=$6

# Unpack Julia
tar -xzf julia-1.9.3-linux-x86_64.tar.gz
tar -xzf netmerge-proj.tar.gz
export PATH=$PWD/julia-1.9.3/bin:$PATH
export JULIA_DEPOT_PATH=$PWD/netmerge-proj

# Run network merging
echo "julia -t${nthreads} --project=netmerge-proj ./perfect_subsets.jl ${netid} ${replicate} ${maxsubsetsize} ${dmethod} ${nsim}
julia -t${nthreads} --project=netmerge-proj ./perfect_subsets.jl ${netid} ${replicate} ${maxsubsetsize} ${dmethod} ${nsim}