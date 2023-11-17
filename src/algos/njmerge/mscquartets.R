library(MSCquartets)
library(ape)

getHBpTable <- function(treefile, cutoff=0.05) {
    gtrees=read.tree(file=treefile)
    taxanames=taxonNames(gtrees)
    QT=quartetTable(gtrees,taxanames)
    RQT=quartetTableResolved(QT)
    pTable=quartetTreeTestInd(RQT,"T3")
    HolmBonferroni(pTable,"T3",cutoff)
}
write.HBpTable <- function(treefile, outfile, cutoff=0.05) {
    write.csv(getHBpTable(treefile, cutoff=cutoff), outfile)
}
args = commandArgs(trailingOnly=T)
write.HBpTable(args[1], args[2])


# write.HBpTable("/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/100-ex2.treefile", "/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/100-ex2.hbptab", cutoff=0.01)

# write.HBpTable("/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/10000-constraints1-major.treefile", "/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/10000-constraints1-major.hbptab", cutoff=0.01)



# tab <- getHBpTable("/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/10000-constraints1.treefile", cutoff=0.01)

# ## Large n80h8 example from `src/tests/netnjmerge/mergealgo.jl`
# file = "/Users/nkolbow/repos/network-merging/src/tests/netnjmerge/10000-sims.treefile"
# gtrees=read.tree(file=file)
# taxanames=taxonNames(gtrees)

# QTstart <- Sys.time()
# QT=quartetTable(gtrees,taxanames)
# QTend <- Sys.time()

# RQTstart <- Sys.time()
# RQT=quartetTableResolved(QT)
# RQTend <- Sys.time()

# pTabStart <- Sys.time()
# pTable=quartetTreeTestInd(RQT,"T3")
# pTabEnd <- Sys.time()

# HBstart <- Sys.time()
# HBpTable=HolmBonferroni(pTable,"T3",.05)
# HBend <- Sys.time()

# print(paste0("QT time: ", QTend-QTstart))
# print(paste0("RQT time: ", RQTend-RQTstart))
# print(paste0("pTable time: ", pTabEnd-pTabStart))
# print(paste0("HBpTable time: ", HBend-HBstart))

# pcutoff=0.05
# nsignif=sum(HBpTable$HBp_T3 <= pcutoff)
# signifQuartets=matrix(NA, nrow=nsignif, ncol=4)
# matrowidx = 1
# i = 1
# while(i <= nsignif && matrowidx <= nsignif) {
#     if(HBpTable$HBp_T3[matrowidx] <= pcutoff) {
#         j = 1
#         for(name in taxanames) {
#             if(HBpTable[i, name] == 1) {
#                 signifQuartets[matrowidx, j] <- name
#                 j <- j + 1
#                 if(j == 5) break
#             }
#         }
#         matrowidx <- matrowidx + 1
#     }
#     i <- i + 1
# }
