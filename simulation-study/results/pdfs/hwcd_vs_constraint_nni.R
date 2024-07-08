library(tidyverse)
library(ggplot2)
source("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/fig-helpers.R")


# Something showing sum of NNIs vs. HWCD(N, hat(N))
p1 <- df %>%
        filter(estRFerror != -1 & netid == "n50r2") %>%
    add_error_bins_zero() %>%
    mutate(est_diff = estRFerror - constraint_error_sum,
           max_subset_size = ordered(max_subset_size)) %>%
    group_by(max_subset_size, gauss_error_level, netid) %>%
    summarise(diff_se = sd(est_diff), diff_mean = mean(est_diff)) %>%
    ggplot(aes(x = max_subset_size, y = diff_mean, fill = gauss_error_level)) +
        geom_col(position = "dodge", alpha = 0.7, color = "black") +
        geom_errorbar(aes(ymin = diff_mean - diff_se, ymax = diff_mean + diff_se), position = "dodge") +
        labs(
            x = "Maximum Constraint Size",
            y = expression(paste("HWCD(N, ", widehat(N), ") - Constraint Errors")),
            fill = "Matrix Error",
            title = "A (N50R2)"
        )

p2 <- df %>%
    filter(estRFerror != -1 & netid == "n100r10") %>%
    add_error_bins_zero() %>%
    mutate(est_diff = estRFerror - constraint_error_sum,
           max_subset_size = ordered(max_subset_size)) %>%
    group_by(max_subset_size, gauss_error_level, netid) %>%
    summarise(diff_se = sd(est_diff), diff_mean = mean(est_diff)) %>%
    ggplot(aes(x = max_subset_size, y = diff_mean, fill = gauss_error_level)) +
        geom_col(position = "dodge", alpha = 0.7, color = "black") +
        geom_errorbar(aes(ymin = diff_mean - diff_se, ymax = diff_mean + diff_se), position = "dodge") +
        labs(
            x = "Maximum Constraint Size",
            y = expression(paste("HWCD(N, ", widehat(N), ") - Constraint Errors")),
            fill = "Matrix Error",
            title = "B (N100R10)"
        )

p3 <- df %>%
    filter(estRFerror != -1 & netid == "n200r20") %>%
    add_error_bins_zero() %>%
    mutate(est_diff = estRFerror - constraint_error_sum,
           max_subset_size = ordered(max_subset_size)) %>%
    group_by(max_subset_size, gauss_error_level, netid) %>%
    summarise(diff_se = sd(est_diff), diff_mean = mean(est_diff)) %>%
    ggplot(aes(x = max_subset_size, y = diff_mean, fill = gauss_error_level)) +
        geom_col(position = "dodge", alpha = 0.7, color = "black") +
        geom_errorbar(aes(ymin = diff_mean - diff_se, ymax = diff_mean + diff_se), position = "dodge") +
        labs(
            x = "Maximum Constraint Size",
            y = expression(paste("HWCD(N, ", widehat(N), ") - Constraint Errors")),
            fill = "Matrix Error",
            title = "C (N200R20)"
        )

p4 <- df %>%
    filter(estRFerror != -1 & netid == "n500r25") %>%
    add_error_bins_zero() %>%
    mutate(est_diff = estRFerror - constraint_error_sum,
           max_subset_size = ordered(max_subset_size)) %>%
    group_by(max_subset_size, gauss_error_level, netid) %>%
    summarise(diff_se = sd(est_diff), diff_mean = mean(est_diff)) %>%
    ggplot(aes(x = max_subset_size, y = diff_mean, fill = gauss_error_level)) +
        geom_col(position = "dodge", alpha = 0.7, color = "black") +
        geom_errorbar(aes(ymin = diff_mean - diff_se, ymax = diff_mean + diff_se), position = "dodge") +
        labs(
            x = "Maximum Constraint Size",
            y = expression(paste("HWCD(N, ", widehat(N), ") - Constraint Errors")),
            fill = "Matrix Error",
            title = "D (N500R25)"
        )

(p1 + p2) / (p3 + p4)
patched <- (p1 + p2) / (p3 + p4)

pdf(file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/hwcd_vs_constraint_nni.pdf", patched, width=12, height=8)
patched
dev.off()
