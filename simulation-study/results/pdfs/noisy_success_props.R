library(tidyverse)
library(ggplot2)
source("/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/fig-helpers.R")



p1 <- plot_success_prop_stacked_bar("n50r2") & labs(title="", x="") & guides(fill="none")
p2 <- plot_success_prop_stacked_bar("n100r5") & labs(title="", x="", y="") & guides(fill="none")
p3 <- plot_success_prop_stacked_bar("n100r10") & labs(title="") & guides(fill="none")
p4 <- plot_success_prop_stacked_bar("n200r10") & labs(title="", y="")

patched <- ((p1 + p2) / (p3 + p4)) + 
    plot_layout(tag_level = "new") +
    plot_annotation(tag_levels = list(c("A (N50R2)", "B (N100R5)", "C (N100R10)", "D (N200R10)")))
patched







pdf(file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/results/noisy_success_props.pdf", width=10, height=9)
patched
dev.off()
