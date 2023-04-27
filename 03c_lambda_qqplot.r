# Load packages------
load_pkg <- rlang::quos(tidyverse, data.table, ggplot2, ggpubr, ggrepel, RColorBrewer, car, speedglm, foreach, doParallel, kableExtra, meta, qqman, patchwork)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

# Load BOLT, meta_wod, meta_wd------
bolt <- fread("/user/work/ac16225/malaria/results/bolt/bolt_neutrophil.african.bgen")
plink_bim <- fread("/user/work/ac16225/malaria/data/genotypes/bed/data.chr1-22.bim")
bolt_filtered <- bolt[SNP %in% plink_bim$V2, ]
meta_wod <- fread("/user/work/ac16225/malaria/results/snptest/res_wod_meta/meta_chr1-22_all.txt")
meta_wd <- fread("/user/work/ac16225/malaria/results/snptest/res_wd_meta/meta_chr1-22_all.txt")
meta_wod_unfiltered <- fread("/user/work/ac16225/malaria/results/snptest/res_wod_meta/meta_chr1-22_all.txt.unfiltered")
meta_wd_unfiltered <- fread("/user/work/ac16225/malaria/results/snptest/res_wd_meta/meta_chr1-22_all.txt.unfiltered")

# Calculate genetic control (lambda)
bolt[, chisq := qchisq(1 - bolt$P_BOLT_LMM_INF, 1)]
ginfl_bolt <- signif( median(bolt$chisq)/qchisq(0.5,1), d = 4 )

bolt_filtered[, chisq := qchisq(1 - bolt_filtered$P_BOLT_LMM_INF, 1)]
ginfl_bolt_f <- signif( median(bolt_filtered$chisq)/qchisq(0.5,1), d = 4 )

meta_wod[, chisq := qchisq(1 - meta_wod$P_value, 1)]
ginfl_meta_wod <- signif( median(meta_wod$chisq)/qchisq(0.5,1), d = 4 )

meta_wd[, chisq := qchisq(1 - meta_wd$P_value, 1)]
ginfl_meta_wd <- signif( median(meta_wd$chisq)/qchisq(0.5,1), d = 4 )

meta_wod_unfiltered[, chisq := qchisq(1 - meta_wod_unfiltered$P_value, 1)]
ginfl_meta_wod_unfiltered <- signif( median(meta_wod_unfiltered$chisq)/qchisq(0.5,1), d = 4 )

meta_wd_unfiltered[, chisq := qchisq(1 - meta_wd_unfiltered$P_value, 1)]
ginfl_meta_wd_unfiltered <- signif( median(meta_wd_unfiltered$chisq)/qchisq(0.5,1), d = 4 )

ginfl_dt <- data.table(
    gwas = c("bolt", "bolt_filtered", "meta_wod", "meta_wd", "meta_wod_unfiltered", "meta_wd_unfiltered"),
    gc_lambda = c(ginfl_bolt, ginfl_bolt_f, ginfl_meta_wod, ginfl_meta_wd, ginfl_meta_wod_unfiltered, ginfl_meta_wd_unfiltered)
)

fwrite(ginfl_dt, "/user/work/ac16225/malaria/results/tables/5-genomic_inflation.txt", quote=F, row.names=F, sep="\t")


# Create variables needed to plot qqplot
# Start with bolt
ci <- 0.95
nSNPs <- nrow(bolt)

bolt[ , `:=` (observed = -log10(sort(bolt$P_BOLT_LMM_INF)), 
    expected = -log10(ppoints(nSNPs)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))))
]
bolt_sub <- bolt %>%
  filter(expected <= 2) %>%
  sample_frac(0.01)

bolt_sup <- bolt %>%
  filter(expected > 2)

bolt_small <- rbind(bolt_sub, bolt_sup)

# meta_wod
ci <- 0.95
nSNPs <- nrow(meta_wod)

meta_wod[ , `:=` (observed = -log10(sort(meta_wod$P_value)), 
    expected = -log10(ppoints(nSNPs)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))))
]
meta_wod_sub <- meta_wod %>%
  filter(expected <= 2) %>%
  sample_frac(0.01)

meta_wod_sup <- meta_wod %>%
  filter(expected > 2)

meta_wod_small <- rbind(meta_wod_sub, meta_wod_sup)

# meta_wd
ci <- 0.95
nSNPs <- nrow(meta_wd)

meta_wd[ , `:=` (observed = -log10(sort(meta_wd$P_value)), 
    expected = -log10(ppoints(nSNPs)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))))
]
meta_wd_sub <- meta_wd %>%
  filter(expected <= 2) %>%
  sample_frac(0.01)

meta_wd_sup <- meta_wd %>%
  filter(expected > 2)

meta_wd_small <- rbind(meta_wd_sub, meta_wd_sup)

# Plot qqplots------
# https://danielroelfs.com/blog/how-i-make-qq-plots-using-ggplot/
bolt_qq <- ggplot(bolt_small, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = "black", size = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
       y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
  theme_pubr()

meta_wod_qq <- ggplot(meta_wod_small, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = "black", size = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
       y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
  theme_pubr()

meta_wd_qq <- ggplot(meta_wd_small, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = "black", size = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
       y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
  theme_pubr()

qqplots <- ggarrange(bolt_qq, meta_wod_qq, meta_wd_qq, common.legend = TRUE, ncol = 1, nrow = 3, labels="AUTO")

ggsave(plot = bolt_qq,
       filename = "/user/work/ac16225/malaria/results/figures/XX_qqplot_main.tiff", 
        device="tiff", width = 270, height = 110, dpi=300, compression = "lzw", units="mm")

ggsave(plot = qqplots,
       filename = "/user/work/ac16225/malaria/results/figures/XX_qqplot_sup.tiff", 
        device="tiff", width = 210, height = 210, dpi=300, compression = "lzw", units="mm")

