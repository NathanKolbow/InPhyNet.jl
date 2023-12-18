#!/bin/bash

# Extract julia & project files
tar -xzf julia-1.9.4-linux-x86_64.tar.gz
tar -xzf netmerge.tar.gz

# add Julia binary to PATH
export PATH=$_CONDOR_SCRATCH_DIR/julia-1.9.4/bin:$PATH
# add Julia packages to DEPOT variable
export JULIA_DEPOT_PATH=$_CONDOR_SCRATCH_DIR/netmerge

julia --project=netmerge --threads=8 pipeline-1i_iv.jl "n40h4" "pipeline-1i_iv.csv" 10000
julia --project=netmerge --threads=8 pipeline-1i_iv.jl "n80h8" "pipeline-1i_iv.csv" 10000
julia --project=netmerge --threads=8 pipeline-1i_iv.jl "n160h16" "pipeline-1i_iv.csv" 10000
julia --project=netmerge --threads=8 pipeline-1i_iv.jl "n240h24-unbalanced" "pipeline-1i_iv.csv" 10000