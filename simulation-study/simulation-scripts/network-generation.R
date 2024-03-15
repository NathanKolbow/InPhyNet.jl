# Run w/ R v4.2.0 (Vigorous Calisthenics)
# Should be run from within directory `network-merging/simulation-study/`
library(ape)
library(SiPhyNetwork)
library(tidyverse)

# Multiple workstations...
basedir <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/"
if(!dir.exists(basedir)) basedir <- "C:\\Users\\Nathan\\repos\\network-merging\\"
if(!dir.exists(basedir)) basedir <- "/Users/Nathan/repos/network-merging/"
basedir <- paste0(basedir, "simulation-study/")

# Import helper functions
source(paste0(basedir, "simulation-scripts/helpers/net-gen-helpers.R"))

top_params <- list(
    num_taxa = c(50, 50, 100, 100, 200, 200, 500, 500, 1000, 1000),
    num_retics = c(2, 5, 5, 10, 10, 20, 25, 50, 50, 100)
)

set.seed(42)
subnet_size <- 25
for(idx in 1:length(top_params$num_taxa)) {
    ntaxa <- top_params$num_taxa[idx]
    nretic <- top_params$num_retics[idx]
    
    output_file <- paste0(
        basedir,
        paste0("data/networks/n", ntaxa, "r", nretic, ".netfile")
    )
    if(file.exists(output_file)) { file.remove(output_file, paste0(output_file, "_copy")) }

    basemsg <- paste0("\rSimulating n", ntaxa, "r", nretic, " networks ")
    cat(paste0(basemsg, "(0/100)"))

    nu <- 0.1
    for(replicate in 1:100) {
        cat(paste0(basemsg, "(", replicate, "/100)"))

        # 1. Generate the tree that the subnetworks will be placed in
        ntips_tre <- ntaxa / subnet_size
        tre <- generate_tree(ntips_tre)

        # 2. Generate the subnetworks & rename them
        subnet_data <- generate_subnets(ntips_tre, subnet_size, nretic, nu)
        subnets <- rename_subnet_tips(subnet_data$subnets, subnet_size)
        nu <- subnet_data$nu

        # 3. Place the subnetworks
        final_newick <- write.net(tre)
        for(idx in 1:ntips_tre) {
            final_newick <- str_replace(final_newick, paste0(tre$tip.label[idx], ":"), paste0("placeholder", idx, ":"))
        }

        retic_idx <- 1
        for(idx in 1:ntips_tre) {
            subnet_newick <- write.net(subnets[[idx]])
            subnet_newick <- str_sub(subnet_newick, 1, str_length(subnet_newick)-1)

            subnet_nretic <- nrow(subnets[[idx]]$reticulation)
            if(subnet_nretic > 0) {
                subnet_newick <- rename_retics(subnet_newick, retic_idx)
                retic_idx <- retic_idx + subnet_nretic
            }

            final_newick <- str_replace(final_newick, paste0("placeholder", idx, ":"), subnet_newick)
        }

        # 4. Save results
        write(final_newick, file=output_file, append=TRUE)
    }
    cat("\n")
}
