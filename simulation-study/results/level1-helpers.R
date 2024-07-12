library(ggplot2)
library(tidyverse)
library(patchwork)
library(cowplot)
library(hash)
library(ggpubr)
options(dplyr.summarise.inform = FALSE)
theme_set(theme_half_open())

GRAD_N3_PALETTE <- c("#deebf7", "#9ecae1", "#51a0d8")
GRAD_N7_PALETTE <- rainbow(7)

# Data
print("Reading D_std0_l1.csv")
df_std0 <- read.csv("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/data/D_std0_l1.csv")
std0_hash <- hash()
for(i in 1:nrow(df_std0)) {
    row <- df_std0[i,]
    std0_hash[[paste0(100 * row$ntaxa + row$rep)]] <- row$std0
}

add_std0 <- function(df) {
    df$std0 <- -1
    ntaxa <- df[1,"numtaxa"]
    for(repnum in 1:100) {
        iter_filt <- df$replicate_num == repnum
        if(any(iter_filt)) {
            df[iter_filt, "std0"] <- std0_hash[[paste0(100 * ntaxa + repnum)]]
        }
    }
    df
}

l1_net_dir <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/data/networks-level1/"
df_level1 <- NULL
for(ntaxa in c(50, 100, 200, 500, 1000, 2500)) {
    print(paste0("Reading n", ntaxa, "-l1.csv"))
    df_level1 <- rbind(
        df_level1,
        read.csv(paste0(l1_net_dir, "n", ntaxa, "-l1.csv")) %>% add_std0()
    )
}
df_level1$gauss_error <- df_level1$std0 / (df_level1$gauss_error + df_level1$std0)


# Figures
plot_hwcd_heatmap_quad_l1 <- function(ntaxa, retic_maxes, m, width = 5, height = 4) {
    ymax <- filter(df_level1, numtaxa == ntaxa & max_subset_size == m) %>%
        select(estRFerror) %>% max()

    plots <- list()
    for(i in 1:4) {
        r_max <- retic_maxes[i] + 1
        r_min <- 1
        if(i > 1) r_min <- retic_maxes[i-1]

        gg_df <- df_level1 %>%
            filter(numtaxa == ntaxa & nretics_true < r_max & nretics_true >= r_min & max_subset_size == m) %>%
            mutate(
                gauss_error_rounded = round(width * gauss_error, digits = 1) / width,
                constraint_error_sum = (round(constraint_error_sum / height + 1e-9, digits=0)) * height,
                estRFerror = ifelse(estRFerror >= 0, estRFerror, NA)
            ) %>%
            group_by(gauss_error_rounded, constraint_error_sum) %>%
            summarise(
                mean_estRFerror = mean(estRFerror, na.rm = TRUE),
                tile_alpha = mean(!is.na(estRFerror))
            ) %>%
            filter(tile_alpha > 0)
        
        tile_width <- sort(unique(gg_df$gauss_error_rounded))[2] - sort(unique(gg_df$gauss_error_rounded))[1]
        gg_ttl <- paste0(LETTERS[i], " (", r_min, " <= r < ", r_max, ")")

        plots[[i]] <- gg_df %>%
            ggplot(aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror, alpha = tile_alpha)) +
                geom_tile(width = tile_width, height = height) +
                scale_fill_gradientn(colors = GRAD_N7_PALETTE, na.value="transparent", limits = c(0, ymax)) +
                scale_alpha_continuous(range = c(0.1, 1), guide = 'none') +
                labs(x = "Distance Matrix Signal", y = "Total Constraint Error (HWCD)", title = gg_ttl, fill = "HWCD")
    }
    
    # top <- ggarrange(plots[[1]], plots[[2]], nrow=1, ncol=2, common.legend = TRUE, legend = 'right')
    # bot <- ggarrange(plots[[3]], plots[[4]], nrow=1, ncol=2, common.legend = TRUE, legend = 'right')
    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow=2, ncol=2, common.legend = TRUE, legend = 'right')
    # ggarrange(top, bot, nrow=2, ncol=1, common.legend = TRUE)
}


plot_hwcd_heatmap_single_l1 <- function(ntaxa, m, r = FALSE, width = 5, height = 4) {
    gg_df <- df_level1 %>%
        filter(numtaxa == ntaxa & max_subset_size == m)
    
    if(!is.logical(r) && is.vector(r)) gg_df <- filter(gg_df, nretics_true %in% r)

    gg_df <- gg_df %>%
        mutate(
            gauss_error_rounded = round(width * gauss_error, digits = 1) / width,
            constraint_error_sum = (round(constraint_error_sum / height + 1e-9, digits=0)) * height,
            estRFerror = ifelse(estRFerror >= 0, estRFerror, NA)
        ) %>%
        group_by(gauss_error_rounded, constraint_error_sum, max_subset_size) %>%
        summarise(
            mean_estRFerror = mean(estRFerror, na.rm = TRUE),
            tile_alpha = mean(!is.na(estRFerror))
        ) %>%
        filter(tile_alpha > 0)
    
    ymax <- gg_df %>% select(mean_estRFerror) %>% max()
    tile_width <- sort(unique(gg_df$gauss_error_rounded))[2] - sort(unique(gg_df$gauss_error_rounded))[1]

    gg_df %>%
        ggplot(aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror, alpha = tile_alpha)) +
            geom_tile(width = tile_width, height = height) +
            scale_fill_gradientn(colors = GRAD_N7_PALETTE, na.value="transparent", limits = c(0, ymax)) +
            scale_alpha_continuous(range = c(0.1, 1)) +
            labs(
                x = "Distance Matrix Signal",
                y = "Total Constraint Error (HWCD)",
                fill = "HWCD",
                title = paste0("Inferred Network Error, n = ", ntaxa),
                alpha = "Success Rate"
            )
}


plot_hwcd_heatmap_m_l1 <- function(ntaxa, width = 5, height = 4, ms = c(5, 10, 15, 20, 25, 30)) {
    temp_df <- df_level1 %>%
        filter(numtaxa == ntaxa)
    ymax <- quantile((temp_df %>% filter(estRFerror >= 0))$estRFerror, 0.999)
    xmin <- quantile((temp_df %>% filter(estRFerror >= 0))$gauss_error, 0.0001)

    plots <- list()
    for(i in 1:length(ms)) {
        gg_df <- temp_df %>%
            filter(max_subset_size == ms[i]) %>%
            mutate(
                gauss_error_rounded = round(width * gauss_error, digits = 1) / width,
                constraint_error_sum = (round(constraint_error_sum / height + 1e-9, digits=0)) * height,
                estRFerror = ifelse(estRFerror >= 0, estRFerror, NA)
            ) %>%
            group_by(gauss_error_rounded, constraint_error_sum) %>%
            summarise(
                mean_estRFerror = mean(estRFerror, na.rm = TRUE),
                tile_alpha = mean(!is.na(estRFerror))
            ) %>%
            filter(tile_alpha > 0)
        
        tile_width <- sort(unique(gg_df$gauss_error_rounded))[2] - sort(unique(gg_df$gauss_error_rounded))[1]
        p <- ggplot(gg_df, aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror, alpha = tile_alpha)) +
                geom_tile(width = tile_width, height = height) +
                scale_fill_gradientn(colors = GRAD_N7_PALETTE, na.value="transparent", limits = c(0, ymax)) +
                scale_alpha_continuous(range = c(0.1, 1)) +
                labs(
                    x = "Distance Matrix Signal",
                    y = "Total Constraint Error (HWCD)",
                    title = paste0(LETTERS[i], " (m = ", ms[i], ")"),
                    fill = "HWCD",
                    alpha = "Success Rate"
                ) +
                scale_y_continuous(limits = c(0, ymax)) +
                scale_x_continuous(limits = c(xmin, 1))
        
        if(i == 1 || i == 2 || i == 4 || i == 5)
            p <- p & guides(fill = 'none', alpha = 'none')
        if(i == 2 || i == 3 || i == 5 || i == 6)
            p <- p & labs(y = NULL)
        
        if(length(ms) == 6) {
            if(i <= 3)
                p <- p & labs(x = NULL)
        }
        
        plots[[i]] <- p
    }
    
    if(length(ms) == 6) {
        return(ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=2, ncol=3, common.legend = TRUE, legend = 'right'))
    } else if(length(ms) == 3) {
        return(ggarrange(plots[[1]], plots[[2]], plots[[3]], nrow=1, ncol=3, common.legend = TRUE, legend = 'right'))
    } else {
        return(plots)
    }
}


plot_est_retics_l1 <- function(ntaxa) {
    gg_df <- df_level1 %>%
        filter(estRFerror >= 0 & numtaxa == ntaxa) %>%
        group_by(nretics_true, max_subset_size) %>%
        summarise(avg_est_retics = mean(nretics_est))
    ymax <- max(gg_df$nretics_true)

    ggplot(gg_df, aes(x = nretics_true, y = avg_est_retics, color = as.factor(max_subset_size))) +
        geom_smooth(se = FALSE, method = "loess", formula = y ~ x, alpha = 1) +
        # geom_line(linewidth = 1) +
        geom_abline(slope = 1, intercept = 0, color = "black", alpha = 0.5, linewidth = 1) +
        labs(
            x = "Number of Reticulations",
            y = "Estimated Reticulations",
            color = "Max Subset Size"
        ) +
        scale_y_continuous(limits = c(0, ymax)) +
        scale_x_continuous(limits = c(0, ymax))
}


plot_hwcd_4nets_l1 <- function(ntaxas, m, ymaxes) {
    length(ntaxas) == 4 || error("Need 4 ntaxa arguments")
    
    heights <- c(4, 4, 10, 12)
    widths <- c(8, 12, 10, 16)
    ps <- list()
    for(i in 1:length(ntaxas)) {
        gg_df <- df_level1 %>%
            filter(numtaxa == ntaxas[i] & max_subset_size == m) %>%
            mutate(
                gauss_error_rounded = round(widths[i] * gauss_error, digits = 1) / widths[i],
                constraint_error_sum = (round(constraint_error_sum / heights[i] + 1e-9, digits=0)) * heights[i],
                estRFerror = ifelse(estRFerror >= 0, estRFerror, NA)
            ) %>%
            group_by(gauss_error_rounded, constraint_error_sum) %>%
            summarise(
                mean_estRFerror = mean(estRFerror, na.rm = TRUE),
                tile_alpha = mean(!is.na(estRFerror))
            ) %>%
            filter(tile_alpha > 0)
        
        ttl <- paste0(LETTERS[i], " (n = ", ntaxas[i], ")")
        tile_width <- sort(unique(gg_df$gauss_error_rounded))[2] - sort(unique(gg_df$gauss_error_rounded))[1]
        p <- gg_df %>%
            ggplot(aes(x = gauss_error_rounded, y = constraint_error_sum, fill = mean_estRFerror, alpha = tile_alpha)) +
                geom_tile(width = 1 / (10 * widths[i]), height = heights[i]) +
                scale_fill_gradientn(colors = GRAD_N7_PALETTE, na.value="transparent", limits = c(0, ymaxes[i])) +
                scale_alpha_continuous(range = c(0.1, 1), guide = 'none') +
                labs(x = "Distance Matrix Signal", y = "Total Constraint Error (HWCD)", title = ttl, fill = "HWCD") +
                scale_y_continuous(limits = c(0, ymaxes[i])) +
                scale_x_continuous(limits = c(0.5, 1.0))
        
        if(i == 2 || i == 4) p <- p & labs(y = "")
        else p <- p & guides(fill = 'none')
        if(i == 1 || i == 2) p <- p & labs(x = "")

        ps[[i]] <- p
    }

    top <- ggarrange(ps[[1]], ps[[2]], nrow=1, ncol=2, common.legend = TRUE, legend = 'right')
    bot <- ggarrange(ps[[3]], ps[[4]], nrow=1, ncol=2, common.legend = TRUE, legend = 'right')
    ggarrange(top, bot, nrow=2, ncol=1, common.legend = FALSE)
}

