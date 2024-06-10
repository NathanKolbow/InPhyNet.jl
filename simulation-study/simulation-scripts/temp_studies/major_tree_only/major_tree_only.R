# RESULT: no statistically significant difference between only major tree & not only major tree

library(ggplot2)
library(tidyverse)


add_error_bins <- function(df) {
    factor_lvls <- c("low", "med", "high", "very high")
    gauss_cuts <- quantile(df$gauss_sd, c(0.1, 0.4, 0.7))
    nni_cuts <- quantile(df$sum_constraint_hwcd, c(0.1, 0.4, 0.7))
    if(!("gauss_error_level" %in% colnames(df))) {
        df <- df %>%
            mutate(
                gauss_error_level = ifelse(gauss_sd < gauss_cuts[1], "low", ifelse(gauss_sd < gauss_cuts[2], "med", ifelse(gauss_sd < gauss_cuts[3], "high", "very high"))),
                gauss_error_level = ordered(gauss_error_level, levels = factor_lvls),
                nni_error_level = ifelse(sum_constraint_hwcd < nni_cuts[1], "low", ifelse(sum_constraint_hwcd < nni_cuts[2], "med", ifelse(sum_constraint_hwcd < nni_cuts[3], "high", "very high"))),
                nni_error_level = ordered(nni_error_level, levels = factor_lvls)
            )
    }
    df
}
df <- read.csv("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/temp_studies/major_tree_only/major_tree_only.csv") %>%
    mutate(major_tree_only = as.logical(major_tree_only))

# Prop success
mean(filter(df, major_tree_only == TRUE)$hwcd == -1)    # 0.548
mean(filter(df, major_tree_only == FALSE)$hwcd == -1)   # 0.542
# after all unrooted --> 0.467, 0.474

# HWCD
n50_df <- df %>%
    filter(hwcd != -1) %>%
    filter(net_id == "n50r2") %>%
    add_error_bins()
n100_df <- df %>%
    filter(hwcd != -1) %>%
    filter(net_id == "n100r10") %>%
    add_error_bins()
n200_df <- df %>%
    filter(hwcd != -1) %>%
    filter(net_id == "n200r10") %>%
    add_error_bins()

ggplot(n50_df, aes(x = gauss_error_level, y = hwcd, color = major_tree_only)) +
    geom_boxplot() +
    facet_grid(m~nni_error_level)
ggplot(n100_df, aes(x = gauss_error_level, y = hwcd, color = major_tree_only)) +
    geom_boxplot() +
    facet_grid(m~nni_error_level)
ggplot(n200_df, aes(x = gauss_error_level, y = hwcd, color = major_tree_only)) +
    geom_boxplot() +
    facet_grid(m~nni_error_level)

ggplot(n50_df, aes(x = sum_constraint_hwcd, y = hwcd, linetype = major_tree_only, color = as.factor(m))) +
    geom_smooth(se = FALSE) +
    geom_abline(slope = 1, intercept = 0, color = "black")

# Major tree only is about the same.
fit_df <- df %>%
    filter(hwcd != -1) %>%
    mutate(major_tree_only = as.numeric(major_tree_only),
           hwcd_diff = hwcd - sum_constraint_hwcd)

fit.n50 <- lm(hwcd_diff ~ gauss_sd + m + major_tree_only, data=filter(fit_df, net_id == "n50r2"))
summary(fit.n50)

fit.n100 <- lm(hwcd_diff ~ gauss_sd + m + major_tree_only, data=filter(fit_df, net_id == "n100r10"))
summary(fit.n100)

fit.n200 <- lm(hwcd_diff ~ gauss_sd + m + major_tree_only, data=filter(fit_df, net_id == "n200r10"))
summary(fit.n200)

summary(
    lm(hwcd ~ net_id * gauss_sd * sum_constraint_hwcd * m + major_tree_only, data=fit_df)
)

# Major tree does [not] improve completion rate
fit_df.n50 <- df %>%
    filter(net_id == "n50r2") %>%
    mutate(success = as.factor(hwcd != -1))

fit.glm.n50 <- glm(success ~ as.factor(m) + major_tree_only, family=binomial(), data=fit_df.n50)
summary(fit.glm.n50)

fit_df.n100 <- df %>%
    filter(net_id == "n100r10") %>%
    mutate(success = hwcd != -1)

fit.glm.n100 <- glm(success ~ as.factor(m) + major_tree_only, family=binomial(), data=fit_df.n100)
summary(fit.glm.n100)

fit_df.n200 <- df %>%
    filter(net_id == "n200r10") %>%
    mutate(success = hwcd != -1)

fit.glm.n200 <- glm(success ~ as.factor(m) + major_tree_only, family=binomial(), data=fit_df.n200)
summary(fit.glm.n200)
