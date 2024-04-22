library(tidyverse)
library(ggplot2)
source("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/fig-helpers.R")


df <- read.df()

##########################################
# FIGURES FOR DATA W/ PERFECT INPUT DATA #
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
plot_success_rate_vs_binned_errors(filter(full_df, netid == "n50r2"))

###################################################################
# HWCD & HWCD W/O MISSING RETIC FIGURES FOR DATA W/ PERTURBATIONS #
###################################################################


# Heatmap for `netid`
plot_hwcd_heatmap(net_df(df, "n50r2"), plot_factor = 2)
plot_hwcd_heatmap(net_df(df, "n50r2"), plot_factor = 2, without_extra_retics = TRUE)

# HWCD w/ all retics vs. only identified retics
plot_compare_hwcd_with_and_wo_identified_retics(net_df(df, "n50r2"))
