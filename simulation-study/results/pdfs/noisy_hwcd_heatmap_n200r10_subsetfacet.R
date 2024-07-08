library(tidyverse)
library(ggplot2)
source("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/fig-helpers.R")


ymax <- max(net_df("n200r10")$estRFerror)

ps <- list()
ms <- c(5, 10, 15, 20, 25, 30)
factors <- c(8, 8, 8, 8,  8, 8)
heights <- c(4, 6, 4, 4, 4, 4)
for(i in 1:6) {
    p <- plot_hwcd_heatmap_std0("n200r10", only_subset_size=ms[i],
        plot_factor=factors[i], without_extra_retics=FALSE, subset_facet=FALSE, col_lims=c(0, ymax), tile_height=heights[i]) +
        scale_y_continuous(limits=c(0, max(net_df("n200r10")$constraint_error_sum))) &
        labs(title=paste0(LETTERS[i], " (m = ", ms[i], ")"))
    if(i < 4) {
        p <- p & labs(x=NULL)
    }
    if(i != 1 && i != 4) {
        p <- p & labs(y=NULL)
    }
    ps[[i]] <- p
}

arr <- ggarrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ps[[6]], nrow=2, ncol=3, common.legend = TRUE, legend = "right")
arr



pdf(file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/noisy_hwcd_heatmap_n200r10_subsetfacet.pdf", width=14, height=9)
arr
dev.off()

# 
