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
    print(paste0("Found nu: ", nu))
    return(nu)
}

findnets <- function(ntaxa, nu, goalhybs) {
    idx <- 1
    nets <- sim.bdh.taxa.ssa(   # placeholder nets
            n=ntaxa, numbsim=1, nu=nu, hybprops=c(0.5, 0, 0.5),
            mu=0, hyb.inher.fxn=make.beta.draw(10, 10),
            lambda=1
        )

    while(idx <= 100) {
        net <- sim.bdh.taxa.ssa(
            n=ntaxa, numbsim=1, nu=nu, hybprops=c(0.5, 0, 0.5),
            mu=0, hyb.inher.fxn=make.beta.draw(10, 10),
            lambda=1
        )
        if(!is.numeric(net[[1]]) && getNetworkLevel(net[[1]]) == 1 && getnhybs(net) == goalhybs) {
            nets[[idx]] <- net[[1]]
            print(idx)
            idx <- idx + 1
        }
    }
    return(nets)
}

# Find networks
findandsavenets <- function(ntaxa, prod) {
    goalhybs <- floor(ntaxa * prod)
    nu <- findnu(ntaxa, goalhybs)
    nets <- findnets(ntaxa, nu, goalhybs)
    write.net(nets, file=paste0("data/networks/n", ntaxa, "h", goalhybs, ".netfile"))
}

# 50 taxa
findandsavenets(50, 0.05)   # done
findandsavenets(50, 0.10)   # done

# 100 taxa
findandsavenets(100, 0.05)
findandsavenets(100, 0.10)

# 200 taxa
findandsavenets(200, 0.05)
findandsavenets(200, 0.10)

# 500 taxa
findandsavenets(500, 0.05)
findandsavenets(500, 0.10)

# 1000 taxa
findandsavenets(1000, 0.05)
findandsavenets(1000, 0.10)