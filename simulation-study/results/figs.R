library(tidyverse)
library(ggplot2)
source("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/fig-helpers.R")


# `df` and `df_std0` are loaded in `fig-helpers.R` and
# are used automatically instead of needing to be passed into each function
#
# df <- read_df()


##########################################
# FIGURES FOR SIMS W/ PERFECT INPUT DATA #
##########################################

no_added_error_df <- filter(full_df, constraint_error_sum == 0 & gauss_error == 0 & estRFerror != -1) %>%
    filter(max_subset_size %in% c(10, 15, 20, 25))  # just so there are fewer things to look @ in the plot pane
ggplot(no_added_error_df, aes(x = netid, y = estRFerror, group = netid)) +
    geom_boxplot() +
    facet_wrap(~max_subset_size)


##################################################
# SUCCESS RATE FIGURES FOR DATA W/ PERTURBATIONS #
##################################################

# x = Gauss error, y = NNI error, color = probability to fail
plot_success_rate_vs_binned_errors(net_df("n50r2"))
plot_success_rate_vs_binned_errors(net_df("n100r10"))
plot_success_rate_vs_binned_errors(net_df("n200r10"))

###################################################################
# HWCD & HWCD W/O MISSING RETIC FIGURES FOR DATA W/ PERTURBATIONS #
###################################################################

# Heatmap for `netid`
plot_hwcd_heatmap("n50r2", plot_factor = 2)
plot_hwcd_heatmap("n50r2", plot_factor = 2, without_extra_retics = TRUE)
plot_hwcd_heatmap("n50r2", plot_factor = 2, without_extra_retics = TRUE, subset_facet = TRUE)

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
plot_hwcd_grouped_boxplots("n200r10", without_extra_retics = TRUE)
plot_hwcd_grouped_boxplots("n500r25", without_extra_retics = TRUE)
plot_hwcd_grouped_boxplots("n1000r100", without_extra_retics = TRUE)

# Best replicate plots
plot_best_replicate_hwcd_grouped_boxplots("n50r2", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n100r5", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n200r10", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n500r25", without_extra_retics = TRUE)
plot_best_replicate_hwcd_grouped_boxplots("n1000r100", without_extra_retics = TRUE)
