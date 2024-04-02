library(tidyverse)
library(ggplot2)

setwd("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/temp_studies")

# Data for `max subset size = 25` and `replicatenum = 1`
df <- read.csv("nni_distribution_data.csv") %>%
    filter(est_errors >= 0) %>%
    mutate(`Total Moves` = as.factor(num_total_moves),
           prop_nets_with_moves = as.factor(paste0(round(100 * num_nets_with_moves / ifelse(net_id == "n100r5", 5, 21), digits=0), "% have NNI")))

# n100r5
pdf("nni_distribution_figs_n100r5.pdf", 12, 5)
df %>% 
    filter(net_id == "n100r5") %>%
    ggplot(aes(x = constraint_diffs, y = est_errors, color = `Total Moves`)) +
    facet_wrap(~ prop_nets_with_moves) +
    geom_jitter(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    labs(x = "Sum of NNI errors", y = "Merged network error") +
    ggtitle("Net: n100r5")
dev.off()

# n500r25
pdf("nni_distribution_figs_n500r25.pdf", 12, 8)
df %>% 
    filter(net_id == "n500r25") %>%
    ggplot(aes(x = constraint_diffs, y = est_errors, color = `Total Moves`)) +
    facet_wrap(~ prop_nets_with_moves) +
    geom_jitter(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    labs(x = "Sum of NNI errors", y = "Merged network error") +
    ggtitle("Net: n500r25")
dev.off()


# Data for `max subset size = 20` and `replicatenum = 2`
df <- read.csv("nni_distribution_data_2.csv") %>%
    filter(est_errors >= 0) %>%
    mutate(`Total Moves` = as.factor(num_total_moves),
           prop_nets_with_moves = as.factor(paste0(round(100 * num_nets_with_moves / ifelse(net_id == "n100r5", 5, 21), digits=0), "% have NNI")))