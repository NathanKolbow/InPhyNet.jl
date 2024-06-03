library(ggplot2)
library(tidyverse)

df <- read.csv("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/data/est_data_output/n200r10.csv")
df$est_newick <- NULL
df$ils_level <- factor(df$ils_level, levels = c("low", "med", "high"))

ggplot(df, aes(x = ils_level, y = est_net_hwcd)) +
    geom_point() +
    geom_point(aes(x = ils_level, y = est_constraint_hwcd_sum), color="red") +
    facet_grid(n_loci~seq_len) +
    labs(x = "ILS", y = "HWCD", title = "top = seq_len, side = ngt")

