library(tidyverse)
library(ggplot2)
source("fig-helpers.R")

df <- read.csv("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/data/est_data_output/n200r10-bad-params.csv") %>%
    mutate(
        est_newick = NULL, distance_method = NULL,
        hwcd = est_net_hwcd, est_net_hwcd = NULL
    )

# Outcomes: hwcd, runtime serial, runtime parallel
# Dependents: n_loci, seq_len, ils_level
ggplot(df, aes(x = n_loci, y = est_net_hwcd_only_identified_retics, color = as.ordered(ils_level), shape = as.ordered(seq_len))) +
    geom_jitter(width=50, height=0, size=8, alpha=0.75, stroke=3) +
    scale_y_continuous(limits=c(0, 300)) +
    facet_wrap(~max_subset_size)

ggplot(df, aes(x = n_loci, y = hwcd, color = as.ordered(ils_level), shape = as.ordered(seq_len))) +
    geom_jitter(width=50, height=0, size=8, alpha=0.75) +
    scale_y_continuous(limits=c(0, 300))

ggplot(df, aes(x = n_loci, y = est_constraint_hwcd_sum, color = as.ordered(ils_level), shape = as.ordered(seq_len))) +
    geom_jitter(width=50, height=0, size=8, alpha=0.75) +
    scale_y_continuous(limits=c(0, 300))

ggplot(df, aes(x = n_loci, y = est_net_hwcd_only_identified_retics - est_constraint_hwcd_sum, color = as.ordered(ils_level), shape = as.ordered(seq_len))) +
    geom_jitter(width=50, height=0, size=8, alpha=0.75) +
    scale_y_continuous(limits=c(0, 300))
