library(tidyverse)
library(ggplot2)

options(dplyr.summarise.inform = FALSE)

result_files <- dir("simulation-study/data/output/")
result_files <- result_files[grepl("n[0-9]*r[0-9]*.csv", result_files)]

full_df <- lapply(
    result_files,
    function(x) {
        read.csv(paste0("simulation-study/data/output/", x))
    }) %>%
    do.call("rbind", .) %>%
    mutate(
        netid = as.factor(paste0("n", numtaxa, "r", nretics_true)),
        constraint_sizes = NULL # not using these for now
    )

# -1 proportion plots
plot_neg_one_props <- function(df, step_size = 0.1, std0=1, title = "") {
    fig_df <- c()

    start_center <- min(df$gauss_error) + step_size / 2
    end_center <- max(df$gauss_error) - step_size / 2

    for(x in seq(start_center, end_center, by = step_size)) {
        iter_df <- filter(df, gauss_error >= x - step_size / 2 & gauss_error < x + step_size / 2) # nolint: error.
        if(nrow(iter_df) == 0) { next }

        fig_df <- rbind(
            fig_df,
            group_by(iter_df, netid, max_subset_size) %>%
                summarize(x = x / std0, y = mean(estRFerror >= 0))
        )
    }

    p <- ggplot(fig_df, aes(x = x, y = y, color = factor(max_subset_size))) +
        geom_point() +
        facet_wrap(~ netid) +
        labs(x = "Noise = Normal(x, x)", "std0 multiplier", y = "Proportion of Runs that Succeeded") +
        ggtitle(title)
    p
}

(p <- plot_neg_one_props(full_df))
(p <- plot_neg_one_props(full_df, step_size = 0.025))


# Accuracy vs. everything (n50r2 & n50r5)
gg_df <- filter(full_df, estRFerror >= 0 & (netid %in% c("n50r2", "n50r5"))) %>%
    sample_n(1e4)


ggplot(gg_df, aes(x = gauss_error, y = constraint_error_sum, color = estRFerror)) +
    geom_jitter(width = 0, height = 1, alpha = 0.5) +
    facet_grid(netid ~ max_subset_size) +
    scale_color_gradientn(
        colors = rainbow(7),
        breaks = c(0, 10, 25, 50, 75, 100)
    ) +
    labs(x = "Noise = Normal(x, x)", y = "Sum of constraint errors (HWCD)") +
    ggtitle("Successful runs only")

gg_df %>%
    filter(gauss_error <= 2) %>%
    ggplot(aes(x = gauss_error, y = constraint_error_sum, color = estRFerror)) +
        geom_jitter(width = 0, height = 1, alpha = 0.5) +
        facet_grid(netid ~ max_subset_size) +
        scale_color_gradientn(colors = rainbow(7)) +
    labs(x = "Noise = Normal(x, x)", y = "Sum of constraint errors (HWCD)") +
    ggtitle("Successful runs only")

# Accuracy vs. everything (n200r10 & n200r20)
gg_df <- filter(full_df, estRFerror >= 0 & (netid %in% c("n200r10", "n200r20")))


ggplot(gg_df, aes(x = gauss_error, y = constraint_error_sum, color = estRFerror)) +
    geom_jitter(width = 0, height = 1, alpha = 0.5) +
    facet_grid(netid ~ max_subset_size) +
    scale_color_gradientn(colors = rainbow(7)) +
    labs(x = "Noise = Normal(x, x)", y = "Sum of constraint errors (HWCD)") +
    ggtitle("Successful runs only")
ggplot(gg_df, aes(x = gauss_error, y = constraint_error_sum, color = majortreeRF)) +
    geom_jitter(width = 0, height = 1, alpha = 0.5) +
    facet_grid(netid ~ max_subset_size) +
    scale_color_gradientn(colors = rainbow(7)) +
    labs(x = "Noise = Normal(x, x)", y = "Sum of constraint errors (HWCD)") +
    ggtitle("Successful runs only")

gg_df %>%
    filter(gauss_error <= 2) %>%
    ggplot(aes(x = gauss_error, y = constraint_error_sum, color = estRFerror)) +
        geom_jitter(width = 0, height = 1, alpha = 0.5) +
        facet_grid(netid ~ max_subset_size) +
        scale_color_gradientn(colors = rainbow(7)) +
    labs(x = "Noise = Normal(x, x)", y = "Sum of constraint errors (HWCD)") +
    ggtitle("Successful runs only")
gg_df %>%
    filter(gauss_error <= 2) %>%
    ggplot(aes(x = gauss_error, y = constraint_error_sum, color = majortreeRF)) +
        geom_jitter(width = 0, height = 1, alpha = 0.5) +
        facet_grid(netid ~ max_subset_size) +
        scale_color_gradientn(colors = rainbow(7)) +
    labs(x = "Noise = Normal(x, x)", y = "Sum of constraint errors (HWCD)") +
    ggtitle("Successful runs only")
