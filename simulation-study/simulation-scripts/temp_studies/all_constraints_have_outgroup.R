library(tidyverse)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)
theme_set(theme_linedraw())

df <- read.csv("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/temp_studies/all_constraints_have_outgroups.csv") %>%
    mutate(max_subset_size_label = paste0("max subset size = ", max_subset_size),
        gauss_error_level_og = ordered(gauss_error_level, levels = c("low", "med", "high")),
        gauss_error_level = ordered(paste0("gauss = ", gauss_error_level), levels = c("gauss = low", "gauss = med", "gauss = high")))

# Proportion of successful runs
# gg_df <- df %>%
#     group_by(netid, replicate_num, max_subset_size_label, gauss_error_level) %>%
#     summarise(prop_runs_good_default = mean(net_hwcd_default >= 0),
#         prop_runs_good_outgroup = mean(net_hwcd_outgroup >= 0))
# ggplot(gg_df, aes(x = replicate_num, y = prop_runs_good_default, shape = "Default", color = "Default")) +
#     geom_point(size = 3) +
#     geom_point(aes(x = replicate_num, y = prop_runs_good_outgroup, shape = "Outgroup", color = "Outgroup"), size = 3) +
#     facet_grid(gauss_error_level ~ netid + max_subset_size_label) +
#     labs(x = "Replicate", y = "Proportion of runs succeeded") +
#     scale_y_continuous(limits = c(0, 1))

# Proportion of successful runs w/ groupings
gg_df <- df %>%
    group_by(netid, replicate_num, max_subset_size_label, gauss_error_level, gauss_error_level_og) %>%
    summarise(prop_runs_good_default = mean(net_hwcd_default >= 0),
        prop_runs_good_outgroup = mean(net_hwcd_outgroup >= 0)) %>%
    mutate(prop_better = prop_runs_good_outgroup - prop_runs_good_default,
        result = ifelse(prop_better < 0, "Default", ifelse(prop_better == 0, "Tie", "All have outgroup")))
ggplot(gg_df, aes(x = gauss_error_level_og, y = prop_better, color = result, shape = result)) +
    geom_point(size = 3) +
    facet_grid(max_subset_size_label ~ netid) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
    labs(x = "Gaussian Noise Level",
        y = "(Proportion default runs succeeded) - (Proportion outgroup runs succeeded)",
        color = "Which succeeds more often?") +
    scale_shape_discrete(guide = "none")

# Prop either succeeds
gg_df <- df %>%
    group_by(netid, replicate_num, max_subset_size_label, gauss_error_level, gauss_error_level_og) %>%
    summarise(prop_runs_good_default = mean(net_hwcd_default >= 0),
        prop_runs_good_outgroup = mean(net_hwcd_outgroup >= 0),
        prop_any_succeed = mean(net_hwcd_default >= 0 | net_hwcd_outgroup >= 0))
ggplot(gg_df, aes(x = gauss_error_level_og)) +
    geom_point(aes(y = prop_runs_good_default, shape = "Default", color = "Default")) +
    geom_point(aes(y = prop_runs_good_outgroup, shape = "Outgroup", color = "Outgroup")) +
    geom_point(aes(y = prop_any_succeed, shape = "Either", color = "Either")) +
    facet_grid(max_subset_size_label ~ netid)

# HWCD
gg_df <- filter(df, net_hwcd_default >= 0 & net_hwcd_outgroup >= 0)
ggplot(gg_df) +
    geom_boxplot(aes(group = constraint_diffs_default, x = constraint_diffs_default, y = net_hwcd_default, color = "Default", shape = "Default"), alpha = 0.75) +
    geom_boxplot(aes(group = constraint_diffs_outgroup, x = constraint_diffs_outgroup, y = net_hwcd_outgroup, color = "Outgroup", shape = "Outgroup"), alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0)

# HWCD differences
gg_df <- filter(df, net_hwcd_default >= 0 & net_hwcd_outgroup >= 0) %>%
    group_by(netid, replicate_num, max_subset_size_label, gauss_error_level, gauss_error_level_og) %>%
    summarise(mean_hwcd_default = mean(net_hwcd_default),
        mean_hwcd_outgroup = mean(net_hwcd_outgroup)) %>%
    mutate(mean_hwcd_diff = mean_hwcd_default - mean_hwcd_outgroup,
        better_with_outgroups = mean_hwcd_diff > 0)
ggplot(gg_df, aes(x = gauss_error_level_og, y = mean_hwcd_diff, color = better_with_outgroups)) +
    geom_point(size = 3) +
    facet_grid(max_subset_size_label ~ netid) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
    labs(x = "Gaussian Noise Level",
        y = "HWCD(default) - HWCD(outgroup)",
        color = "Better w/ outgroups?")
