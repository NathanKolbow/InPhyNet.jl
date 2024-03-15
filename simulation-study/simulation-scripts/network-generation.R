# Should be run from within directory `network-merging/simulation-study/`
library(ape)
library(SiPhyNetwork)
source("helpers/net-gen-helpers.R")

# Multiple workstations...
basedir <- "/mnt/ws/home/nkolbow/repos/network-merging/"
if(!dir.exists(basedir)) basedir <- "C:\\Users\\Nathan\\repos\\network-merging\\"
if(!dir.exists(basedir)) basedir <- "/Users/Nathan/repos/network-merging/"


top_params <- list(
    num_taxa = c(50, 50, 100, 100, 200, 200, 500, 500, 1000, 1000),
    num_retics = c(2, 5, 5, 10, 10, 20, 25, 50, 50, 100)
)

for(idx in 1:length(top_params$num_taxa)) {
    ntaxa <- top_params$num_taxa[idx]
    nretic <- top_params$num_retics[idx]

    # 1. Generate the tree that the subnetworks will be placed in


    # 2. Generate the subnetworks


    # 3. Place the subnetworks


    # 4. Save results
    
}