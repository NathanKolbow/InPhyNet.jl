library(tidyverse)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)
theme_set(theme_linedraw())

df <- read.csv("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/temp_studies/all_constraints_have_outgroups.csv") %>%
    mutate(max_subset_size_label = paste0("m = ", max_subset_size),
        gauss_error_level = ordered(paste0("gauss = ", gauss_error_level), levels = c("gauss = low", "gauss = med", "gauss = high")))

# Proportion of successful runs
gg_df <- df %>%
    group_by(netid, replicate_num, max_subset_size_label, gauss_error_level) %>%
    summarise(prop_runs_good_default = mean(net_hwcd_default >= 0),
        prop_runs_good_outgroup = mean(net_hwcd_outgroup >= 0))

ggplot(gg_df, aes(x = replicate_num, y = prop_runs_good_default, shape = "Default", color = "Default")) +
    geom_point(size = 3) +
    geom_point(aes(x = replicate_num, y = prop_runs_good_outgroup, shape = "Outgroup", color = "Outgroup"), size = 3) +
    facet_grid(gauss_error_level ~ netid + max_subset_size_label) +
    labs(x = "Replicate", y = "Proportion of runs succeeded") +
    scale_y_continuous(limits = c(0, 1))

# HWCD
gg_df <- filter(df, net_hwcd_default >= 0 & net_hwcd_outgroup >= 0)
ggplot(gg_df) +
    geom_boxplot(aes(group = constraint_diffs_default, x = constraint_diffs_default, y = net_hwcd_default, color = "Default", shape = "Default"), alpha = 0.75) +
    geom_boxplot(aes(group = constraint_diffs_outgroup, x = constraint_diffs_outgroup, y = net_hwcd_outgroup, color = "Outgroup", shape = "Outgroup"), alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0)

gg_df <- filter(df, net_hwcd_default >= 0 & net_hwcd_outgroup >= 0) %>%
    group_by(netid, replicate_num, max_subset_size_label, gauss_error_level) %>%
    summarise(mean_hwcd_default = mean(net_hwcd_default),
        mean_hwcd_outgroup = mean(net_hwcd_outgroup)) %>%
    mutate(mean_hwcd_diff = mean_hwcd_default - mean_hwcd_outgroup,
        better_with_outgroups = mean_hwcd_diff > 0)
ggplot(gg_df, aes(x = gauss_error_level, y = mean_hwcd_diff, color = better_with_outgroups)) +
    geom_point() +
    facet_grid(netid ~ max_subset_size_label) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
    annotate("text", x = 1, y = 5, label = ">0: better w/ outgroups", color = "red")

