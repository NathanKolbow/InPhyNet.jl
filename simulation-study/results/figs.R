library(tidyverse)
library(ggplot2)

options(dplyr.summarise.inform = FALSE)

read.df <- function() {
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
}
full_df <- read.df()


# Accuracy vs. everything (n50r2 & n50r5)
net_hwcd_heatmap <- function(df, title="Successful runs only", do_sample_n=TRUE, width=0.35, height=2) {
    if(do_sample_n)
        df <- df %>% sample_n(min(1.5e4, nrow(df)))

    ggplot(df, aes(x = gauss_error, y = constraint_error_sum, fill = majortreeRF)) +
        geom_tile(width=width, height=height, alpha=0.45) +
        facet_grid(netid ~ max_subset_size) +
        scale_fill_gradientn(
            colors = rainbow(7),
            # breaks = c(0, 10, 25, 50, 75, 100)
        ) +
        labs(x = "Noise = Normal(x, x)", y = "Sum of constraint errors (HWCD)") +
        ggtitle(title)
}

# Need to change to gauss_error == 0 when that data arrives...
no_added_error_df <- filter(full_df, constraint_error_sum == 0 & gauss_error == 0 & estRFerror != -1) %>%
    filter(max_subset_size %in% c(10, 15, 20, 25))  # just so there are fewer things to look @ in the plot pane
ggplot(no_added_error_df, aes(x = netid, y = estRFerror, group = netid)) +
    geom_boxplot() +
    facet_wrap(~max_subset_size)


# x = Gauss error, y = NNI error, color = probability to fail
factor_lvls <- c("low", "med", "high", "very high")
n200r10_df <- filter(full_df, netid == "n200r10") %>%
    mutate(
        gauss_error_level = ifelse(gauss_error < 0.75, "low", ifelse(gauss_error < 2.5, "med", ifelse(gauss_error < 4.5, "high", "very high"))),
        gauss_error_level = ordered(gauss_error_level, levels = factor_lvls),
        nni_error_level = ifelse(constraint_error_sum < 25, "low", ifelse(constraint_error_sum < 50, "med", ifelse(constraint_error_sum < 65, "high", "very high"))),
        nni_error_level = ordered(nni_error_level, levels = factor_lvls)
    )
gg_df <- n200r10_df %>%
    group_by(gauss_error_level, nni_error_level) %>%
    summarise(mean_val = round(100 * mean(estRFerror != -1), digits=2)) %>%
    mutate(prop = mean_val / 100)

ggplot(gg_df, aes(x = gauss_error_level, y = nni_error_level, fill = mean_val, label = paste0(mean_val, "%"))) +
    geom_tile() +
    geom_text() +
    labs(x = "Distance Matrix Noise Level", y = "Constraint Network Noise Level", fill = "% runs succeeded") +
    scale_fill_gradient(limits = c(0, 100)) +
    ggtitle("Percent runs succeeded based on input noise")


# Heatmap for `netid`
gg_df <- filter(full_df, netid == "n500r25" & estRFerror != -1) %>%
    mutate(gauss_error_rounded = round(4 * gauss_error, digits = 1) / 4) %>%
    group_by(gauss_error_rounded, constraint_error_sum) %>%
    summarise(mean_estRFerror = mean(estRFerror))

ggplot(gg_df, aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror)) +
    geom_tile(width = 0.1, height = 2) +
    scale_fill_gradientn(limits = c(0, max(gg_df$mean_estRFerror)), colors = rainbow(7))
