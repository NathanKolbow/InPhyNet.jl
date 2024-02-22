# Should be in directory `network-merging/simulation-study/`

library(ape)
library(SiPhyNetwork)

# Multiple workstations...
basedir <- "/mnt/ws/home/nkolbow/repos/network-merging/"
if(!dir.exists(basedir)) basedir <- "C:\\Users\\Nathan\\repos\\network-merging\\"
if(!dir.exists(basedir)) basedir <- "/Users/Nathan/repos/network-merging/"

# Helper functions
getnhybs <- function(ssa_nets) {
    nhyblist <- c()
    for(i in 1:length(ssa_nets)) {
        if(!is.numeric(ssa_nets[[i]])) {
            nhybs <- dim(ssa_nets[[i]]$reticulation)[1]
            nhyblist <- c(nhyblist, nhybs)
        }
    }
    return(nhyblist)
}

findnu <- function(ntaxa, goalhybs) {
    nu <- 0.1
    ssa_nets_hybs <- c(goalhybs + 100)
    while(!(mean(ssa_nets_hybs) < goalhybs + 1)) {
        nu <- nu - nu / 10
        ssa_nets <- sim.bdh.taxa.ssa(
            n=ntaxa, numbsim=100, nu=nu, hybprops=c(0.5, 0, 0.5),
            mu=0, hyb.inher.fxn=make.beta.draw(10, 10),
            lambda=1
        )
        ssa_nets_hybs <- getnhybs(ssa_nets)
    }
    return(nu)
}

findnets <- function(ntaxa, nu, goalhybs, maxlevel) {
    idx <- 1
    nets <- sim.bdh.taxa.ssa(   # placeholder nets
            n = ntaxa, numbsim = 1, nu = nu, hybprops = c(0.5, 0, 0.5),
            mu = 0, hyb.inher.fxn = make.beta.draw(10, 10),
            lambda = 1
        )

    while (idx <= 100) {
        net <- sim.bdh.taxa.ssa(
            n = ntaxa, numbsim = 1, nu = nu, hybprops = c(0.5, 0, 0.5),
            mu = 0, hyb.inher.fxn = make.beta.draw(10, 10),
            lambda = 1
        )

        if(!is.numeric(net[[1]]) && getNetworkLevel(net[[1]]) <= maxlevel && getnhybs(net) == goalhybs) {
            nets[[idx]] <- net[[1]]
            cat(paste0("\rFound ", idx, " networks."))
            idx <- idx + 1
        }
    }
    cat("\n")
    return(nets)
}

telap <- function(st_ti) round(difftime(Sys.time(), st_ti, units = "secs"), 2)
findlowestlevel <- function(ntaxa, nu, goalhybs, searchiters = 1e4, shouldprint = T) {
    minlvl <- Inf
    minlvlcount <- 0
    start_time <- Sys.time()
    total_generated_nets <- 0

    while (total_generated_nets < searchiters) {
        net <- sim.bdh.taxa.ssa(
            n = ntaxa, numbsim = 1, nu = nu, hybprops = c(0.5, 0, 0.5),
            mu = 0, hyb.inher.fxn = make.beta.draw(10, 10),
            lambda = 1
        )
        total_generated_nets <- total_generated_nets + 1

        if(!is.numeric(net[[1]]) && getnhybs(net) == goalhybs) {
            currlvl <- getNetworkLevel(net[[1]])
            if(currlvl < minlvl) {
                minlvl <- currlvl
                minlvlcount <- 1
            } else if(currlvl == minlvl) {
                minlvlcount <- minlvlcount + 1
            }
        }

        if(shouldprint && total_generated_nets %% 100 == 0) {
            t_elap <- round(telap(start_time), 2)
            cat(paste0("\r", total_generated_nets, " generated in ", t_elap, " seconds "))
            est_total_time <- round(t_elap / (total_generated_nets / searchiters), 2)
            cat(paste0("(est. ", est_total_time, "s == ", round(est_total_time / 60, 2), "m)                "))
        }
    }
    if(shouldprint) {
        cat(paste0("\nSearch ", total_generated_nets, " in ", round(telap(start_time) / 60, 2), " minutes.\n"))
    }
    return(minlvl)
}

# Find networks
findandsavenets <- function(ntaxa, prod, maxlevel = 1, searchiters = -1) {
    set.seed(42)

    goalhybs <- floor(ntaxa * prod)
    cat("\rSearching for nu value...")
    nu <- findnu(ntaxa, goalhybs)
    cat(paste0("\rFound nu value: ", nu, "                        \n"))

    if(maxlevel == 1 && searchiters != -1) {
        cat(paste0("Searching for appropriate maximum level with ", searchiters, " iterations...\n"))
        maxlevel <- findlowestlevel(ntaxa, nu, goalhybs, searchiters = searchiters)
        print(paste0("Setting max level: ", maxlevel))
    }

    nets <- findnets(ntaxa, nu, goalhybs, maxlevel = maxlevel)
    write.net(nets, file=paste0(basedir, "simulation-study/data/networks/n", ntaxa, "r", goalhybs, ".netfile"))
}

# !!!!!!!!!!!!!!!!!!!!!
# !! REPRODUCIBILITY !!
# !!!!!!!!!!!!!!!!!!!!!
#
# To reproduce study results, this file should not be
# run as an Rscript "as-is". First, run the code about,
# then run a *single line* of the code below. I.e., if
# you want to generate the n500r50 networks, run the
# code above, then *immediately* run:
# `findandsavenets(500, 0.10, searchiters = 5e3)`


# 50 taxa
findandsavenets(50, 0.05)   # can get level-1 easily
findandsavenets(50, 0.10)   # can get level-1 easily

# 100 taxa
findandsavenets(100, 0.05)  # can get level-1 relatively easily
findandsavenets(100, 0.10, searchiters = 1e5)

# 200 taxa
findandsavenets(200, 0.05, searchiters = 1e4)
findandsavenets(200, 0.10, searchiters = 1e4)

# 500 taxa
findandsavenets(500, 0.05, searchiters = 5e3)
findandsavenets(500, 0.10, searchiters = 5e3)

# 1000 taxa
findandsavenets(1000, 0.05, searchiters = 1e3)
findandsavenets(1000, 0.10, searchiters = 1e3)