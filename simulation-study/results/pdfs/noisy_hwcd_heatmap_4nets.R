library(tidyverse)
library(ggplot2)
library(ggpubr)
source("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/fig-helpers.R")


# HWCD heatmap
m <- 25
max.p1p2 <- net_df("n100r5") %>%
    rbind(net_df("n100r10")) %>%
    filter(max_subset_size == m) %>%
    select(estRFerror) %>%
    max()
max.p3p4 <- net_df("n200r10") %>%
    rbind(net_df("n500r25")) %>%
    filter(max_subset_size == m) %>%
    select(estRFerror) %>%
    max()

max.p1p2 <- (round(max.p1p2 / 25, digits=0) + 1) * 25
max.p3p4 <- (round(max.p3p4 / 25, digits=0) + 1) * 25

p1 <- plot_hwcd_heatmap_std0("n100r5", only_subset_size=m, plot_factor = 8, without_extra_retics = FALSE, subset_facet = FALSE, col_lims=c(0, max.p1p2)) &
    labs(title="A (n = 100, r = 5)", x="")

p2 <- plot_hwcd_heatmap_std0("n100r10", only_subset_size=m, plot_factor = 7, tile_height=3, without_extra_retics = FALSE, subset_facet = FALSE, col_lims=c(0, max.p1p2)) &
    labs(title="B (n = 100, r = 10)", y="", x="")

p3 <- plot_hwcd_heatmap_std0("n200r10", only_subset_size=m, plot_factor = 8, tile_height=4, without_extra_retics = FALSE, subset_facet = FALSE, col_lims=c(0, max.p3p4)) &
    labs(title="C (n = 200, r = 10)")

p4 <- plot_hwcd_heatmap_std0("n500r25", only_subset_size=m, plot_factor = 10, tile_height = 7, without_extra_retics = FALSE, subset_facet = FALSE, col_lims=c(0, max.p3p4)) &
    labs(title="D (n = 500, r = 25)", y="")


arr1 <- ggarrange(p1, p2, nrow=1, ncol=2, common.legend = TRUE, legend = "right")
arr2 <- ggarrange(p3, p4, nrow=1, ncol=2, common.legend = TRUE, legend = "right")
arr <- ggarrange(arr1, arr2, nrow=2, ncol=1)
arr



pdf(file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/noisy_hwcd_heatmap_4nets.pdf", width=14, height=9)
arr
dev.off()
