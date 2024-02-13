# Should be in directory `network-merging/simulation-study/`
library(ape)
library(SiPhyNetwork)

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

findlowestlevel <- function(ntaxa, nu, goalhybs, searchtimeseconds = 120, shouldprint = F) {
    minlvl <- Inf
    minlvlcount <- 0
    start_time <- Sys.time()
    telap <- function(st_ti) round(difftime(Sys.time(), st_ti, units = "secs"), 2)

    while (telap(start_time) < searchtimeseconds) {
        net <- sim.bdh.taxa.ssa(
            n = ntaxa, numbsim = 1, nu = nu, hybprops = c(0.5, 0, 0.5),
            mu = 0, hyb.inher.fxn = make.beta.draw(10, 10),
            lambda = 1
        )

        if(!is.numeric(net[[1]]) && getnhybs(net) == goalhybs) {
            currlvl <- getNetworkLevel(net[[1]])
            if(currlvl < minlvl) {
                minlvl <- currlvl
                minlvlcount <- 1
                if(shouldprint) {
                    cat(paste0("\rNew min level: ", minlvl,
                            " (", minlvlcount, ")",
                            " (", telap(start_time), "s taken)       "))
                }
            } else if(currlvl == minlvl && shouldprint) {
                minlvlcount <- minlvlcount + 1
                cat(paste0("\rNew min level: ", minlvl,
                        " (", minlvlcount, ")",
                        " (", telap(start_time), "s taken)       "))
            }
        }
    }
    if(shouldprint) cat("\n")
    return(minlvl)
}
# findlowestlevel(200, 0.00202755595904453, 20, shouldprint = T)

# Find networks
findandsavenets <- function(ntaxa, prod, maxlevel = 1, searchtimeseconds = -1) {
    goalhybs <- floor(ntaxa * prod)
    cat("\rSearching for nu value...")
    nu <- findnu(ntaxa, goalhybs)
    cat(paste0("\rFound nu value: ", nu, "                        \n"))

    if(maxlevel == 1 && searchtimeseconds != -1) {
        cat(paste0("\rSearching for appropriate maximum level in ", searchtimeseconds, "s..."))
        maxlevel <- findlowestlevel(ntaxa, nu, goalhybs, searchtimeseconds = searchtimeseconds)
        cat(paste0("\rSetting max level: ", maxlevel,
                   "                           \n"))
    }

    nets <- findnets(ntaxa, nu, goalhybs, maxlevel = maxlevel)
    write.net(nets, file=paste0("/mnt/ws/home/nkolbow/repos/network-merging/simulation-study/data/networks/n", ntaxa, "r", goalhybs, ".netfile"))
}

# 50 taxa
findandsavenets(50, 0.05)                   # done
findandsavenets(50, 0.10)                   # done

# 100 taxa
findandsavenets(100, 0.05)                  # done
findandsavenets(100, 0.10, maxlevel = 6)    # done

# 200 taxa
findandsavenets(200, 0.05, maxlevel = 6)    # done
findandsavenets(200, 0.10, searchtimeseconds = 240)

# 500 taxa
findandsavenets(500, 0.05, searchtimeseconds = 240)
findandsavenets(500, 0.10, searchtimeseconds = 240)

# 1000 taxa
findandsavenets(1000, 0.05, searchtimeseconds = 240)
findandsavenets(1000, 0.10, searchtimeseconds = 240)