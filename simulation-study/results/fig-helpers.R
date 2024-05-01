library(ggplot2)
library(tidyverse)
library(hash)
library(patchwork)
library(cowplot)
options(dplyr.summarise.inform = FALSE)
theme_set(theme_half_open())

# Color palettes
GRAD_N3_PALETTE <- c("#deebf7", "#9ecae1", "#51a0d8")
GRAD_N7_PALETTE <- rainbow(7)
#c('#b2182b','#ef8a62','#fddbc7','#f7f7f7','#d1e5f0','#67a9cf','#2166ac')

########################
##### Data loading #####
########################

base_dir <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/"

# Read all the csv files into a single df
net_df <- function(netid_filter) { filter(df, netid == netid_filter) }
net_std0 <- function(netid_filter, repnum) { filter(df_std0, net_id == netid_filter & replicate_num == repnum)$D_std[1] }
net_no_noise_df <- function(netid_filter) { filter(no_noise_df, netid == netid_filter) }

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
            constraint_sizes = NULL,    # not using these for now
            all_have_outgroup = as.logical(all_have_outgroup),
            outgroup_removed_after_reroot = as.logical(outgroup_removed_after_reroot)
        ) %>%
        filter()
    
    # Stitch the std0 values in
    # print("Stitching std0")
    # time_start <- Sys.time()
    # full_df$std0 <- -1
    # std0_hash <- hash()
    # for(i in 1:nrow(df_std0)) {
    #     row <- df_std0[i,]
    #     std0_hash[[paste0(row$net_id, "_", row$replicate_num)]] <- row$D_std
    # }

    # for(netid in unique(full_df$netid)) {
    #     idx_filt <- full_df$netid == netid
    #     for(repnum in 1:100) {
    #         iter_filt <- idx_filt & full_df$replicate_num == repnum
    #         if(any(iter_filt)) {
    #             full_df[iter_filt, ]$std0 <- std0_hash[[paste0(netid, "_", repnum)]]
    #         }
    #     }
    # }
    # elaps <- difftime(Sys.time(), time_start, units="secs")

    print("Skipping stitch step")
    full_df
}
df <- read_df()
no_noise_df <- filter(df, gauss_error == 0 & constraint_error_sum == 0)


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

plot_hwcd_grouped_boxplots <- function(netid, without_extra_retics = FALSE, subset_facet = FALSE, best_replicate = FALSE, best_n_replicates = Inf) {
    gg_df <- net_df(netid) %>%
        filter(estRFerror != -1) %>%
        add_error_bins()

    fig_title <- netid
    if(best_replicate) {
        best_rep <- -1
        best_rep_median <- Inf
        for(rep in unique(gg_df$replicate_num)) {
            rep_med <- filter(gg_df, replicate_num == rep & gauss_error_level == "very high" & nni_error_level == "low")
            rep_med <- ifelse(without_extra_retics, rep_med$esterror_without_missing_retics, rep_med$estRFerror)

            if(rep_med < best_rep_median) {
                best_rep_median <- rep_med
                best_rep <- rep
            }
        }
        gg_df <- filter(gg_df, replicate_num == best_rep)
        fig_title <- paste0(netid, " best replicate (", best_rep, ")")
    } else if (best_n_replicates < Inf) {
        uq_reps <- unique(gg_df$replicate_num)
        max_reps <- min(length(uq_reps), best_n_replicates)

        best_n_reps <- c()
        best_n_reps_medians <- c()
        for(rep in uq_reps) {
            rep_med <- filter(gg_df, replicate_num == rep & gauss_error_level == "very high" & nni_error_level == "low")
            rep_med <- ifelse(without_extra_retics, rep_med$esterror_without_missing_retics, rep_med$estRFerror)
            
            if(length(best_n_reps) < max_reps) {
                best_n_reps <- c(best_n_reps, rep)
                best_n_reps_medians <- c(best_n_reps_medians, rep_med)
            } else {
                less_than <- rep_med < best_n_reps_medians
                if(any(less_than)) {
                    replace_idx <- which(best_n_reps_medians == max(best_n_reps_medians))
                    best_n_reps[replace_idx] <- rep
                    best_n_reps_medians[replace_idx] <- rep_med
                }
            }
        }
        gg_df <- filter(gg_df, replicate_num %in% best_n_reps)
        fig_title <- paste0(netid, " best ", max_reps, "/", length(uq_reps), " replicates")
    }

    if(without_extra_retics) {gg_df$yval <- gg_df$esterror_without_missing_retics}
    else {gg_df$yval <- gg_df$estRFerror}
    ylabel <- ifelse(without_extra_retics, "HWCD(truth w/o missing retics, est)", "HWCD(truth, est)")
    
    ymax <- ifelse(netid == "n100r5", 200,
            ifelse(netid == "n100r10", 200,
            ifelse(netid == "n200r10", 325,
            max(gg_df$yval))))

    ggplot(gg_df, aes(x = gauss_error_level, 
                      y = yval,
                      fill = nni_error_level)) +
        geom_boxplot() +
        scale_y_continuous(limits = c(0, ymax)) +
        labs(x = "Gaussian Noise Level", y = ylabel, title = fig_title)
}


plot_best_replicate_hwcd_grouped_boxplots <- function(netid, ...) {
    plot_hwcd_grouped_boxplots(netid, best_replicate = TRUE, ...)
}


plot_best_p_replicates_hwcd_grouped_boxplots <- function(netid, p, ...) {
    temp_df <- net_df(netid) %>%
        filter(estRFerror != -1) %>%
        add_error_bins()
    total_reps <- length(unique(temp_df$replicate_num))
    choose_n <- p * total_reps
    choose_n <- ifelse(choose_n %% 1 == 0, choose_n, choose_n - (choose_n %% 1) + 1)
    
    plot_hwcd_grouped_boxplots(netid, best_n_replicates = max(1, choose_n), ...)
}


plot_est_hwcd_vs_sum_input_error_line <- function(netid, without_extra_retics = FALSE, ...) {
    gg_df <- filter(net_df(netid), estRFerror != -1) %>%
        add_error_bins() %>%
        sample_n(50000)
    if(without_extra_retics) gg_df$est_error <- gg_df$esterror_without_missing_retics
    else gg_df$est_error <- gg_df$estRFerror
    gg_df <- gg_df[order(gg_df$constraint_error_sum),]

    ggplot(gg_df, aes(x = constraint_error_sum, y = est_error, color = gauss_error_level)) +
        geom_smooth(se = FALSE) +
        geom_abline(slope = 1, intercept = 0) +
        scale_y_continuous(limits = c(0, max(gg_df$est_error)))
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


plot_hwcd_heatmap <- function(netid, plot_factor = 2, tile_width = 1 / (plot_factor * 10), tile_height = 2, without_extra_retics = FALSE, subset_facet = FALSE, use_std0 = FALSE, mark_minimal = FALSE, draw_failures = FALSE, only_subset_size = NA, ...) {
    minimal_error_cutoff <- 4

    gg_df <- net_df(netid)
    if(!draw_failures) {
        gg_df <- filter(gg_df, estRFerror != -1)
    }
    if(!is.na(only_subset_size)) {
        gg_df <- filter(gg_df, max_subset_size == only_subset_size)
    }
    gg_df$estRFerror[gg_df$estRFerror < 0] <- NA
    gg_df$esterror_without_missing_retics[gg_df$esterror_without_missing_retics < 0] <- NA
    gg_df$estRFerror_compare <- gg_df$estRFerror
    gg_df$estRFerror_compare[gg_df$estRFerror_compare < 0] <- Inf

    if(use_std0) gg_df$gauss_error <- gg_df$gauss_error / gg_df$std0
    gg_df <- gg_df %>%
        mutate(gauss_error_rounded = round(plot_factor * gauss_error, digits = 1) / plot_factor)
    
    # Group by max subset size as well if we want to facet it
    if(subset_facet) {
        gg_df <- gg_df %>%
            group_by(gauss_error_rounded, constraint_error_sum, max_subset_size)
    } else {
        gg_df <- gg_df %>%
            group_by(gauss_error_rounded, constraint_error_sum)
    }
    gg_df <- gg_df %>%
            summarise(
                mean_estRFerror = mean(estRFerror, na.rm = TRUE),
                mean_estRFerror_wo_retics = mean(esterror_without_missing_retics, na.rm = TRUE),
                any_with_minimal_error = any(estRFerror_compare <= minimal_error_cutoff), 1, 0)
            #      %>%
            # mutate(minimal_error_alpha = ifelse(any_with_minimal_error, 1, 0))

    p <- NULL
    # Plot either HWCD w/ all retics or w/o unidentified retics
    if(without_extra_retics) {
        p <- ggplot(gg_df, aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror_wo_retics)) +
            geom_tile(width = tile_width, height = tile_height, ...) +
            scale_fill_gradientn(limits = c(0, max(gg_df$mean_estRFerror)),
            colors = GRAD_N7_PALETTE, na.value = "lightgray")
    } else {
        p <- ggplot(gg_df, aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror)) +
            geom_tile(width = tile_width, height = tile_height, ...) +
            scale_fill_gradientn(limits = c(0, max(gg_df$mean_estRFerror)),
            colors = GRAD_N7_PALETTE, na.value = "lightgray")
    }

    if(mark_minimal) {
        p <- p +
            geom_point(aes(alpha = minimal_error_alpha), color = "black", shape = 4) +
            guides(alpha = "none")
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


plot_success_rate_vs_binned_errors <- function(netid) {
    gg_df <- net_df(netid) %>%
        add_error_bins() %>%
        group_by(gauss_error_level, nni_error_level) %>%
        summarise(mean_val = round(100 * mean(estRFerror != -1), digits=2)) %>%
        mutate(prop = mean_val / 100)

    ggplot(gg_df, aes(x = gauss_error_level, y = nni_error_level, fill = mean_val, label = paste0(mean_val, "%"))) +
        geom_tile() +
        geom_text() +
        labs(x = "Distance Matrix Noise Level", y = "Constraint Network Noise Level", fill = "% runs succeeded") +
        scale_fill_gradientn(limits = c(0, 100), colors = GRAD_N3_PALETTE) +
        ggtitle("Percent runs succeeded based on input noise")
}


compare_outgroup_plot_succ_vs_binned_errors <- function(netid) {
    gg_df <- net_df(netid) %>%
        add_error_bins() %>%
        group_by(gauss_error_level, nni_error_level, all_have_outgroup, outgroup_removed_after_reroot) %>%
        summarise(mean_val = round(100 * mean(estRFerror != -1), digits=2)) %>%
        mutate(prop = mean_val / 100,
               xval = paste0(str_sub(paste0(all_have_outgroup), 1, 1), str_sub(paste0(outgroup_removed_after_reroot), 1, 1)))

    ggplot(gg_df, aes())

    gg_df %>%
        ggplot(aes(x = xval, y = prop, color = xval)) +
        geom_abline(slope = 0, intercept = c(0.0, 0.5, 1.0), linetype = "dashed", alpha = 0.5) +
        geom_bar(stat = "identity", position = "dodge", fill = "white", linewidth = 1.5, alpha = 0.9) +
        labs(x = "Distance Matrix Noise Level", y = "Constraint Network Noise Level", fill = "% runs succeeded") +
        facet_grid(nni_error_level ~ gauss_error_level)
}
compare_outgroup_plot_succ_vs_binned_errors("n100r5")


plot_success_rate_vs_binned_errors_lineplot <- function(netid) {
    gg_df <- net_df(netid) %>%
        add_error_bins() %>%
        group_by(gauss_error_level, nni_error_level) %>%
        summarise(mean_val = round(100 * mean(estRFerror != -1), digits=2)) %>%
        mutate(prop = mean_val / 100)

    ggplot(gg_df, aes(x = gauss_error_level, y = prop, group = nni_error_level, color = nni_error_level)) +
        geom_point() +
        geom_line() +
        labs(x = "Gaussian error level",
            y = "Proportion of runs succeeded",
            color = "Constraint error level")
}


plot_success_prop_stacked_bar <- function(netid, no_noise = FALSE) {
    gg_df <- c()
    if(no_noise) gg_df <- net_no_noise_df(netid)
    else gg_df <- net_df(netid)

    gg_df <- gg_df %>%
        add_error_bins() %>%
        group_by(gauss_error_level, nni_error_level) %>%
        summarise(prob = round(100 * mean(estRFerror != -1), digits=2)) %>%
        mutate(stack_prob = prob)

    for(nni_level in levels(gg_df$nni_error_level)) {
        idx <- gg_df$nni_error_level == nni_level

        # gauss == very high
        curr_idx <- idx & gg_df$gauss_error_level == "very high"
        stack_total <- gg_df$stack_prob[curr_idx]

        # gauss == high
        curr_idx <- idx & gg_df$gauss_error_level == "high"
        gg_df$stack_prob[curr_idx] <- gg_df$prob[curr_idx] - stack_total
        stack_total <- gg_df$prob[curr_idx]

        # gauss == med
        curr_idx <- idx & gg_df$gauss_error_level == "med"
        gg_df$stack_prob[curr_idx] <- gg_df$prob[curr_idx] - stack_total
        stack_total <- gg_df$prob[curr_idx]

        # gauss == low
        curr_idx <- idx & gg_df$gauss_error_level == "low"
        gg_df$stack_prob[curr_idx] <- gg_df$prob[curr_idx] - stack_total
        stack_total <- gg_df$prob[curr_idx]
    }

    ggplot(gg_df, aes(x = nni_error_level, y = stack_prob, fill = gauss_error_level, group = gauss_error_level)) +
        geom_bar(stat = "identity", color = "black") +
        geom_text(aes(label=paste0(prob, "%")), position = position_stack(vjust = 0.5), color = "white") +
        labs(x = "Constraint Error Level", y = "Success Probability (%)",
            fill = "Gaussian Noise Level")
}


plot_no_noise_hwcd <- function(without_extra_retics = FALSE) {
    gg_df <- no_noise_df %>%
        filter(estRFerror != -1) %>%
        group_by(max_subset_size, netid) %>%
        summarise(mean_RF = mean(estRFerror))
    
    ggplot(gg_df, aes(x = max_subset_size, y = mean_RF, color = netid)) +
        geom_line() +
        geom_point() +
        scale_x_continuous(breaks=c(5, 10, 15, 20, 25, 30)) +
        geom_abline(slope=0, intercept=0, linetype="dashed", color="black", alpha=0.5)
}


plot_no_noise_success_props <- function(netid) {
    gg_df <- no_noise_df %>%
        group_by(netid) %>%
        summarise(n_success = sum(estRFerror != -1),
                  prop_success = mean(estRFerror != -1))

    ggplot(gg_df, aes(x = netid, y = prop_success)) +
        geom_bar(stat = "identity", color = "black", fill = "white") +
        scale_y_continuous(limits = c(0, 1.05)) +
        geom_text(aes(label = round(prop_success, digits = 2), y = prop_success + 0.02)) +
        labs(x = "Network ID", y = "Proportion of Replicates Succeeded")
}
