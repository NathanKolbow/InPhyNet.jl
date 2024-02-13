library(MSCquartets)
library(ape)
library(tidyverse)

treefile <- "C:\\Users\\Nathan\\repos\\network-merging\\src\\tests\\netnjmerge\\10000-ex2.treefile"
gtrees <- read.tree(file=treefile)
taxanames=taxonNames(gtrees)

runtime <- function(ntrees, ntaxa) {
    trees <- gtrees[1:ntrees]
    starttime = Sys.time()
    QT <- quartetTable(trees, taxanames[1:ntaxa])
    RQT <- quartetTableResolved(QT)
    pTable <- quartetTreeTestInd(RQT, "T3")
    hbpTable <- HolmBonferroni(pTable, "T3", 0.01)
    Sys.time() - starttime
}

results <- c()
for(ntaxa in 5:17) {
    for(ntrees in c(100, 500, 1000)) {
        timetaken <- runtime(ntrees, ntaxa)
        results <- rbind(results, tibble(ntaxa=ntaxa, ntrees=ntrees, rtime=timetaken))
    }
}
results$rtime <- as.numeric(results$rtime)
write.table(results, "mscqruntime.tab")

library(ggplot2)
ggplot(results, aes(y=rtime, x=ntrees)) +
    geom_point()
ggplot(results, aes(y=rtime, x=ntaxa, color=as.character(ntrees))) +
    geom_point()
