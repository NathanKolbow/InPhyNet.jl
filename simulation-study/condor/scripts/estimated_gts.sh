#!/bin/bash
# $1: netid
# $2: replicate number
# $3: nloci
# $4: seq len
# $5: ILS level
# $6: maximum subset size

net_id=$1
replicate=$2
nloci=$3
seq_len=$4
ils_level=$5
m=$6

# Unpack Julia
tar -xzf julia-1.9.3-linux-x86_64.tar.gz
tar -xzf inphynet-proj.tar.gz
export PATH=$PWD/julia-1.9.3/bin:$PATH
export JULIA_DEPOT_PATH=$PWD/inphynet-proj/

# 1. Simulate true gene trees, convert to sequences, then estimate gene trees
./julia-1.9.3/bin/julia --project=$PWD/inphynet-proj/ --optimize=3 -t2 -p4 ./est-gts.jl ${net_id} ${replicate} ${nloci} ${seq_len} ${ils_level} ${m}

# 2. Infer constraints
./julia-1.9.3/bin/julia --project=$PWD/inphynet-proj/ --optimize=3 -t2 -p4 ./infer-constraints.jl ${net_id} ${replicate} ${nloci} ${seq_len} ${ils_level} ${m}

# 3. InPhyNet inference
./julia-1.9.3/bin/julia --project=$PWD/inphynet-proj/ --optimize=3 -t2 -p4 ./inphynet.jl ${net_id} ${replicate} ${nloci} ${seq_len} ${ils_level} ${m}


# Run the pipeline
echo "./julia-1.9.3/bin/julia -t2 -p4 --project=$PWD/inphynet-proj/ ./estimated_gts.jl ${net_id} ${replicate} ${nloci} ${seq_len} ${ils_level} ${m}"
./julia-1.9.3/bin/julia --project=$PWD/inphynet-proj/ --pkgimages=no --optimize=3 -t2 -p4 ./estimated_gts.jl ${net_id} ${replicate} ${nloci} ${seq_len} ${ils_level} ${m}

