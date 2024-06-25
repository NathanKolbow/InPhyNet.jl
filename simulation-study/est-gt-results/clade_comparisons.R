library(ggplot2)
library(ggtree)
library(tidyverse)
library(awtools)

# 1. Only take major trees and set all edge lengths to -1.
tru <- treeio::read.tree(text="((((((((t174,t175),t155),(((t160,t161),t159),(t170,t171))),((((t152,t153),t158),t151),(((t172,t173),(t154,((t162,t166),(t163,(t164,t165))))),((t157,(t168,t169)),(t156,t167))))),(((t179,(t183,((t199,t200),(t191,t192)))),((t193,t194),((t180,t181),(t182,(((t188,t189),t185),t184))))),(((t177,(t178,(t197,t198))),((t190,(t195,t196)),(t186,t187))),t176))),(((((t11,(t12,t13)),(((t14,t15),(t16,((t22,t23),t21))),t8)),(t6,t7)),(t19,t20)),(t1,((((t4,t5),t3),(((t9,t10),(t17,t18)),(t24,t25))),t2)))),((((((t56,(t72,t73)),(t55,(t64,t65))),((((t54,(t68,t69)),t51),(t60,t61)),(t59,(t74,t75)))),((t52,(t70,t71)),((t62,t63),((t57,t58),(t53,(t66,t67)))))),((((((t92,t93),(t83,t84)),(t78,t79)),t77),(t76,(((t90,t91),((t81,t82),((t99,t100),t94))),((t85,t86),t80)))),(((t88,t89),(t97,t98)),((t95,t96),t87)))),(((((t106,t107),((t102,(((t118,t119),((t108,t109),((t122,t123),(t114,t115)))),t104)),((t105,(((t124,t125),t117),(t110,(t120,t121)))),t116))),((((t112,t113),t111),t103),t101)),(((((t145,t146),t142),t140),(t126,((t132,t133),((((t149,t150),t130),(t138,t139)),t128)))),(((((t134,t135),((t141,(t143,t144)),(t147,t148))),t129),t127),(t131,(t136,t137))))),(((t26,((t45,t46),t34)),t29),((((t43,t44),((t33,t41),t30)),(t37,t38)),(((t28,(t35,((t47,t48),t42))),(t31,((t49,t50),t32))),(t27,((t39,t40),t36)))))))),OUTGROUP);")

est <- treeio::read.tree(text="(OUTGROUP,((((t1,(((((t10,t9),(t17,t18)),(t24,t25)),((t3,t5),t4)),t2)),((((t11,(t12,t13)),(((t14,t15),(t16,((t21,t23),t22))),t8)),(t6,t7)),(t19,t20))),((((t151,((t152,t158),t153)),(((t156,t167),(t157,(t168,t169))),((t154,((t162,(t163,(t164,t165))),t166)),(t172,t173)))),((t155,(t174,t175)),(((t159,t161),t160),(t170,t171)))),((t176,((t177,(t178,(t197,t198))),((t186,t187),(t190,(t195,t196))))),((t179,(t183,((t191,(t199,t200)),t192))),(((t180,t181),(t182,(t184,(t185,(t188,t189))))),(t193,t194)))))),(((((((((t100,(t94,t99)),t81),(t82,(t90,t91))),(t80,(t85,t86))),t76),(t77,(t78,(t79,((t83,t84),(t92,t93)))))),((t87,(t95,t96)),((t88,t89),(t97,t98)))),((((((t54,t51),(t68,t69)),(t60,t61)),(t59,(t74,t75))),((t55,(t64,t65)),(t56,(t72,t73)))),((((t52,t71),t70),(((t53,t67),t66),(t57,t58))),(t62,t63)))),(((((t101,(t103,(t111,t113))),t112),(((((t102,(t104,(((t115,(t108,t109)),(((t114,t119),t118),t123)),t122))),t107),((t105,((t110,(t120,t121)),(t125,t124))),t117)),t106),t116)),(((t126,((t128,((t130,(t149,t150)),(t138,t139))),(t132,t133))),((t140,(t145,t146)),t142)),(t127,(((t129,((t134,t135),(((t141,t144),t143),(t147,t148)))),(t136,t137)),t131)))),(((t26,t29),(t34,(t45,t46))),(((((t30,(t33,t41)),(t43,t44)),(t37,t38)),(t42,(t47,t48))),(t28,((((t27,t36),(t39,t40)),((t31,(t49,t50)),t32)),t35))))))));")

tip_categories <- tibble(
    tip_name = c(paste0("t", 1:200), "OUTGROUP"),
    clade_number = c(((1:200 - 1) %/% 25), "OUTGROUP")
)

plot_colored_tree <- function(tre, title = "", hilight_nodes=c(1), hflip=FALSE) {
    # Find the MRCA node number for each clade
    mrcas <- c()
    for(clade_num in unique(tip_categories$clade_number)) {
        temp_names <- tip_categories$tip_name[tip_categories$clade_number == clade_num]
        if(length(temp_names) == 1) { next } # outgroup

        temp_mrca <- temp_names[1]
        for(i in 2:length(temp_names)) {
            mm <- MRCA(tre, temp_names[1], temp_names[i])
            if(length(mm) == 0) { next }
            temp_mrca <- min(temp_mrca, mm)
        }
        mrcas <- c(mrcas, temp_mrca)
    }
    mrcas <- as.integer(mrcas)
    names(mrcas) <- 1:8
    tre <- groupClade(tre, mrcas)

    # Plot stuff
    colors <- awtools::ppalette
    p <- ggtree(tre, aes(color = group)) %<+% tip_categories +
        geom_tiplab(
            color = "black",
            geom = "label",
            label.padding = unit(0.15, "lines"),
            label.size = 0) +
        scale_color_manual(values = c("black", colors))
    
    for(i in 1:length(mrcas)) {
        p <- p +
            geom_hilight(
                node=as.integer(mrcas[i]), alpha = 0.4,
                color = colors[i], fill = colors[i],
                type = "roundrect"
            )
    }
    p <- p & guides(color = "none")
    
    p
}

plot_colored_tree(tru, title = "True major tree", hflip=TRUE)
plot_colored_tree(est, title = "Estimated major tree")


# 2. Save as PDF
# 3. Add lines in Krita
# 4. Add reticulations in Krita
pdf(file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/est-gt-results/tru_clade_comparison.pdf", width=18, height=40)
plot_colored_tree(tru, title = "True major tree")
dev.off()

pdf(file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/est-gt-results/est_clade_comparison.pdf", width=18, height=40)
plot_colored_tree(est, title = "Estimated major tree")
dev.off()



# BOOK: https://yulab-smu.top/treedata-book/chapter5.html
# Line coloring from this plot may also be preferable: https://yulab-smu.top/treedata-book/chapter13.html#hpv58