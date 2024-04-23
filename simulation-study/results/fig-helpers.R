library(ggplot2)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)

########################
##### Data loading #####
########################

base_dir <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/"

# Read all the csv files into a single df
net_df <- function(netid_filter) { filter(df, netid == netid_filter) }
net_std0 <- function(netid_filter, repnum) { filter(df_std0, net_id == netid_filter & replicate_num == repnum)$D_std[1] }

read_std0_df <- function() {
    read.csv(file.path(base_dir, "data", "D_std0.csv"))
}
df_std0 <- read_std0_df()

read_df <- function() {
    output_dir <- file.path(base_dir, "data", "output")
    result_files <- dir(output_dir)
    result_files <- result_files[grepl("n[0-9]*r[0-9]*.csv", result_files)]
    full_df <- lapply(
        result_files,
        function(x) {
            read.csv(paste0(file.path(output_dir, x)))
        }) %>%
        do.call("rbind", .) %>%
        mutate(
            netid = as.factor(paste0("n", numtaxa, "r", nretics_true)),
            constraint_sizes = NULL # not using these for now
        )
    
    # Stitch the std0 values in
    full_df$std0 <- -1

    print("Stitching std0 into df")
    for(netid in unique(full_df$netid)) {
        for(repnum in 1:100) {
            idx_filt <- full_df$netid == netid & full_df$replicate_num == repnum
            if(any(idx_filt)) {
                full_df[idx_filt, ]$std0 <- net_std0(netid, repnum)
            }
        }
    }
    full_df
}
df <- read_df()


add_error_bins <- function(df) {
    factor_lvls <- c("low", "med", "high", "very high")
    gauss_cuts <- quantile(df$gauss_error, c(0.1, 0.4, 0.7))
    nni_cuts <- quantile(df$constraint_error_sum, c(0.1, 0.4, 0.7))
    if(!("gauss_error_level" %in% colnames(df))) {
        df <- df %>%
            mutate(
                gauss_error_level = ifelse(gauss_error < gauss_cuts[1], "low", ifelse(gauss_error < gauss_cuts[2], "med", ifelse(gauss_error < gauss_cuts[3], "high", "very high"))),
                gauss_error_level = ordered(gauss_error_level, levels = factor_lvls),
                nni_error_level = ifelse(constraint_error_sum < nni_cuts[1], "low", ifelse(constraint_error_sum < nni_cuts[2], "med", ifelse(constraint_error_sum < nni_cuts[3], "high", "very high"))),
                nni_error_level = ordered(nni_error_level, levels = factor_lvls)
            )
    }
    df
}


###################
##### Figures #####
###################

plot_hwcd_grouped_boxplots <- function(netid, without_extra_retics = FALSE, subset_facet = FALSE, best_replicate = FALSE) {
    gg_df <- net_df(netid) %>%
        filter(estRFerror != -1) %>%
        add_error_bins()

    fig_title <- netid
    if(best_replicate) {
        best_rep <- -1
        best_rep_median <- Inf
        for(rep in unique(gg_df$replicate_num)) {
            rep_med <- median(filter(gg_df, replicate_num == rep & gauss_error_level == "very high" & nni_error_level == "low")$esterror_without_missing_retics)

            if(rep_med < best_rep_median) {
                best_rep_median <- rep_med
                best_rep <- rep
            }
        }
        gg_df <- filter(gg_df, replicate_num == best_rep)
        fig_title <- paste0(netid, " best replicate (", best_rep, ")")
    }

    if(without_extra_retics) {gg_df$yval <- gg_df$esterror_without_missing_retics}
    else {gg_df$yval <- gg_df$estRFerror}
    ylabel <- ifelse(without_extra_retics, "HWCD(truth w/o missing retics, est)", "HWCD(truth, est)")
    
    ggplot(gg_df, aes(x = gauss_error_level, 
                      y = yval,
                      fill = nni_error_level)) +
        geom_boxplot() +
        labs(x = "Gaussian Noise Level", y = ylabel, title = fig_title)
}


plot_best_replicate_hwcd_grouped_boxplots <- function(netid, ...) {
    plot_hwcd_grouped_boxplots(netid, best_replicate = TRUE, ...)
}


plot_compare_hwcd_with_and_wo_identified_retics <- function(gg_df, sample_size = 1e3) {
    gg_df %>%
        sample_n(sample_size) %>%
        ggplot(aes(x = estRFerror, y = esterror_without_missing_retics)) +
            geom_jitter(width = 0.5, height = 0.5, alpha = 0.85) +
            geom_abline(slope = 1, intercept = 0) +
            labs(x = "HWCD(truth, estimated)", y = "HWCD(truth w/o unidentified retics, estimated)",
                title = "HWCD(true network, estimated network) with vs. without unidentified retics")
}


plot_hwcd_heatmap <- function(netid, plot_factor = 2, tile_width = 1 / (plot_factor * 10), tile_height = 2, without_extra_retics = FALSE, subset_facet = FALSE, use_std0 = FALSE, ...) {
    gg_df <- filter(net_df(netid), estRFerror != -1)
    if(use_std0) gg_df$gauss_error <- gg_df$gauss_error / gg_df$std0
    gg_df <- gg_df %>%
        mutate(gauss_error_rounded = round(plot_factor * gauss_error, digits = 1) / plot_factor)
    
    # Group by max subset size as well if we want to facet it
    if(subset_facet) {
        gg_df <- gg_df %>%
            group_by(gauss_error_rounded, constraint_error_sum, max_subset_size) %>%
            summarise(mean_estRFerror = mean(estRFerror), mean_estRFerror_wo_retics = mean(esterror_without_missing_retics))
    } else {
        gg_df <- gg_df %>%
            group_by(gauss_error_rounded, constraint_error_sum) %>%
            summarise(mean_estRFerror = mean(estRFerror), mean_estRFerror_wo_retics = mean(esterror_without_missing_retics))
    }

    p <- NULL
    # Plot either HWCD w/ all retics or w/o unidentified retics
    if(without_extra_retics) {
        p <- ggplot(gg_df, aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror_wo_retics)) +
            geom_tile(width = tile_width, height = tile_height, ...) +
            scale_fill_gradientn(limits = c(0, max(gg_df$mean_estRFerror)), colors = rainbow(7))
    } else {
        p <- ggplot(gg_df, aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror)) +
            geom_tile(width = tile_width, height = tile_height, ...) +
            scale_fill_gradientn(limits = c(0, max(gg_df$mean_estRFerror)), colors = rainbow(7))
    }

    # Facet if we want to facet
    if(subset_facet) {
        p <- p + facet_wrap(~ max_subset_size)
    }

    # Add labels
    p + labs(
        x = ifelse(use_std0, "std0 Multiplier", "Gaussian Noise = N(x, x)"),
        y = "Sum of Constraint Errors (HWCD)",
        title = ifelse(without_extra_retics, "HWCD(true net w/o unidentified retics, est net)", "HWCD(true net, est net)"),
        fill = "HWCD"
    )
}


plot_hwcd_heatmap_std0 <- function(netid, ...) {
    plot_hwcd_heatmap(netid, use_std0 = TRUE, ...)
}


plot_success_rate_vs_binned_errors <- function(gg_df) {
    gg_df <- add_error_bins(gg_df)

    gg_df <- gg_df %>%
        group_by(gauss_error_level, nni_error_level) %>%
        summarise(mean_val = round(100 * mean(estRFerror != -1), digits=2)) %>%
        mutate(prop = mean_val / 100)

    p <- ggplot(gg_df, aes(x = gauss_error_level, y = nni_error_level, fill = mean_val, label = paste0(mean_val, "%"))) +
        geom_tile() +
        geom_text() +
        labs(x = "Distance Matrix Noise Level", y = "Constraint Network Noise Level", fill = "% runs succeeded") +
        scale_fill_gradient(limits = c(0, 100)) +
        ggtitle("Percent runs succeeded based on input noise")
    return(p)
}