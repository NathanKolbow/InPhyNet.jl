library(tidyverse)
library(ggplot2)
# source("simulation-study/results/fig-helpers.R")
source("simulation-study/results/level1-helpers.R")

# plot_hwcd_heatmap("n200r10", only_subset_size = 25)

(df_level1 %>% filter(numtaxa == 2500) %>% select(nretics_true))[,1] %>% unique() %>% sort()

# n2500
plot_hwcd_heatmap_single_l1(2500, 25)

# n1000
ntaxa <- 1000
retic_maxes <- c(8, 15, 22, Inf)
plot_hwcd_heatmap_single_l1(ntaxa, 25)
plot_hwcd_heatmap_quad_l1(ntaxa, retic_maxes, 25, height=8, width=10)

# n500
ntaxa <- 500
retic_maxes <- c(5, 10, 15, Inf)
plot_hwcd_heatmap_quad_l1(ntaxa, retic_maxes, 25)
plot_hwcd_heatmap_m_l1(ntaxa)

# n200
ntaxa <- 200
retic_maxes <- c(3, 6, 9, Inf)
plot_hwcd_heatmap_quad_l1(ntaxa, retic_maxes, 25)

# n100
ntaxa <- 100
retic_maxes <- c(2, 4, 6, Inf)
plot_hwcd_heatmap_quad_l1(ntaxa, retic_maxes, 25)

# n50
ntaxa <- 50
retic_maxes <- c(Inf)
plot_hwcd_heatmap_single_l1(ntaxa, 25)


# 4 nets
plot_hwcd_4nets_l1(c(100, 200, 500, 1000), 25, c(125, 125, 400, 400))
    


# How many retics estimated?
plot_est_retics_l1(50)
plot_est_retics_l1(100)
plot_est_retics_l1(200)
plot_est_retics_l1(500)





# SVGs
fig_dir <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/level1-pdfs/"

svg(paste0(fig_dir, "hwcd_m_n500.svg"), width=16, height=6)
plot_hwcd_heatmap_m_l1(500, width = 6, height = 8, ms = c(5, 15, 25))
dev.off()

svg(paste0(fig_dir, "retics_est_n500.svg"), width=12, height=8)
plot_est_retics_l1(500)
dev.off()

svg(paste0(fig_dir, "hwcd_4nets_m25.svg"), width=12, height=8)
plot_hwcd_4nets_l1(c(100, 200, 500, 1000), 25, c(125, 125, 400, 400))
dev.off()

svg(paste0(fig_dir, "hwcd_n1000_retics.svg"), width=12, height=8)
plot_hwcd_heatmap_quad_l1(1000, c(8, 15, 23, 50), 25, height=12, width=10)
dev.off()

svg(paste0(fig_dir, "hwcd_n2500.svg"), width = 8, height = 6)
plot_hwcd_heatmap_single_l1(2500, 25, width = 10, height = 16)
dev.off()

svg(paste0(fig_dir, "hwcd_n1000.svg"), width = 8, height = 6)
plot_hwcd_heatmap_single_l1(1000, 25, width = 10, height = 10)
dev.off()
