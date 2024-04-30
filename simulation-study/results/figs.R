library(tidyverse)
library(ggplot2)
library(cowplot)
source("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/fig-helpers.R")


# `df` and `df_std0` are loaded in `fig-helpers.R` and
# are used automatically instead of needing to be passed into each function
#
# df <- read_df()

####################
# ALWAYS OUTGROUP? #
####################

compare_outgroup_plot_succ_vs_binned_errors("n100r10")

# No major completion success difference...
net_df("n100r10") %>%
    add_error_bins() %>%
    group_by(netid, replicate_num, max_subset_size, gauss_error_level, all_have_outgroup, outgroup_removed_after_reroot) %>%
    summarise(prop = mean(estRFerror >= 0)) %>%
    mutate(xval = factor(paste0(str_sub(paste0(all_have_outgroup), 1, 1), str_sub(paste0(outgroup_removed_after_reroot), 1, 1)))) %>%
    ggplot(aes(x = gauss_error_level, y = prop, color = xval, shape = xval)) +
        geom_jitter(size = 3, height = 0, width = 0.25, alpha = 0.85) +
        facet_grid(max_subset_size ~ netid) +
        geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
        labs(x = "Gaussian Noise Level",
            y = "(Proportion default runs succeeded) - (Proportion outgroup runs succeeded)",
            color = "Which succeeds more often?") +
        scale_shape_discrete(guide = "none")

# No major HWCD differences...
net_df("n100r10") %>%
    filter(estRFerror != -1) %>%
    add_error_bins() %>%
    group_by(netid, replicate_num, max_subset_size, gauss_error_level, all_have_outgroup, outgroup_removed_after_reroot) %>%
    mutate(xval = factor(paste0(str_sub(paste0(all_have_outgroup), 1, 1), str_sub(paste0(outgroup_removed_after_reroot), 1, 1)))) %>%
    ggplot(aes(x = xval, y = estRFerror, color = xval)) +
        geom_boxplot() +
        facet_grid(nni_error_level ~ gauss_error_level) +
        geom_abline(slope = 0, intercept = 0, linetype = "dashed")

# No major differences in mean (or median) HWCDs...
mean(filter(net_df("n100r10"), !all_have_outgroup & estRFerror != -1)$estRFerror)
mean(filter(net_df("n100r10"), all_have_outgroup & outgroup_removed_after_reroot & estRFerror != -1)$estRFerror)
mean(filter(net_df("n100r10"), all_have_outgroup & !outgroup_removed_after_reroot & estRFerror != -1)$estRFerror)

mean(filter(net_df("n100r10"), !all_have_outgroup)$estRFerror == -1)
mean(filter(net_df("n100r10"), all_have_outgroup & outgroup_removed_after_reroot)$estRFerror == -1)
mean(filter(net_df("n100r10"), all_have_outgroup & !outgroup_removed_after_reroot)$estRFerror == -1)

##########################################
# FIGURES FOR SIMS W/ PERFECT INPUT DATA #
##########################################

plot_no_noise_success_props("n50r2")

# HWCD for all nets compared across subset sizes
plot_no_noise_hwcd()

##################################################
# SUCCESS RATE FIGURES FOR DATA W/ PERTURBATIONS #
##################################################

# x = Gauss error, y = NNI error, color = probability to fail
plot_success_rate_vs_binned_errors("n50r2")
plot_success_rate_vs_binned_errors("n100r10")
plot_success_rate_vs_binned_errors("n200r10")

plot_success_rate_vs_binned_errors_lineplot("n50r2")
plot_success_rate_vs_binned_errors_lineplot("n100r5")

# another type of graph
plot_success_prop_stacked_bar("n50r2")
plot_success_prop_stacked_bar("n200r10")

###################################################################
# HWCD & HWCD W/O MISSING RETIC FIGURES FOR DATA W/ PERTURBATIONS #
###################################################################

# Heatmap for `netid`
plot_hwcd_heatmap("n50r2", plot_factor = 2)
plot_hwcd_heatmap("n50r2", plot_factor = 2, without_extra_retics = TRUE)
plot_hwcd_heatmap("n50r2", plot_factor = 2, without_extra_retics = TRUE, subset_facet = TRUE)
plot_hwcd_heatmap("n200r10", plot_factor = 2, without_extra_retics = TRUE, subset_facet = TRUE)
plot_hwcd_heatmap("n200r10", plot_factor = 2, without_extra_retics = TRUE, subset_facet = FALSE)
# plot_hwcd_heatmap("n200r10", plot_factor = 2, without_extra_retics = TRUE, mark_minimal = TRUE, subset_facet = TRUE)

# Heatmap w/ std0 as x-axis
plot_hwcd_heatmap_std0("n50r2", plot_factor = 10, without_extra_retics = TRUE)
plot_hwcd_heatmap_std0("n50r5", plot_factor = 10, without_extra_retics = TRUE)
plot_hwcd_heatmap_std0("n100r5", plot_factor = 10, without_extra_retics = TRUE)
plot_hwcd_heatmap_std0("n100r10", plot_factor = 10, without_extra_retics = TRUE)

# HWCD w/ all retics vs. only identified retics
plot_compare_hwcd_with_and_wo_identified_retics(net_df(df, "n50r2"))

# HWCD vs binned noise levels boxplots
plot_hwcd_grouped_boxplots("n50r2", without_extra_retics = TRUE)
plot_hwcd_grouped_boxplots("n100r5", without_extra_retics = TRUE)
plot_hwcd_grouped_boxplots("n100r10", without_extra_retics = TRUE)
plot_hwcd_grouped_boxplots("n200r10", without_extra_retics = TRUE)
plot_hwcd_grouped_boxplots("n500r25", without_extra_retics = TRUE)
plot_hwcd_grouped_boxplots("n1000r100", without_extra_retics = TRUE)

# Best replicate plots
plot_best_replicate_hwcd_grouped_boxplots("n50r2", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n100r5", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n100r10", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n200r10", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n500r25", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n1000r100", without_extra_retics = TRUE)

# Best (100*p)% of replicates plot
plot_best_p_replicates_hwcd_grouped_boxplots("n100r5", p = 0.5, without_extra_retics = TRUE)
plot_best_p_replicates_hwcd_grouped_boxplots("n100r5", p = 0.25, without_extra_retics = TRUE)
plot_best_p_replicates_hwcd_grouped_boxplots("n100r5", p = 0.10, without_extra_retics = TRUE)
plot_best_p_replicates_hwcd_grouped_boxplots("n100r5", p = 0.0001, without_extra_retics = TRUE)

plot_best_p_replicates_hwcd_grouped_boxplots("n200r10", p = 0.5, without_extra_retics = TRUE)
plot_best_p_replicates_hwcd_grouped_boxplots("n200r10", p = 0.25, without_extra_retics = TRUE)
plot_best_p_replicates_hwcd_grouped_boxplots("n200r10", p = 0.10, without_extra_retics = TRUE)
plot_best_p_replicates_hwcd_grouped_boxplots("n200r10", p = 0.0001, without_extra_retics = TRUE)

# Estimated net error vs. input error, binned by Gaussian noise
plot_est_hwcd_vs_sum_input_error_line("n100r5", without_extra_retics = TRUE)
plot_est_hwcd_vs_sum_input_error_line("n100r10", without_extra_retics = TRUE)

############################
# ESTIMATED DATA HWCD FIGS #
############################

