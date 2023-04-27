#!/usr/bin/env Rscript

# Load relevant libraries------
list_of_packages <- c(
  "data.table",
  "TwoSampleMR",
  "ieugwasr",
  "tidyr",
  "ggplot2", "dplyr", "ggforestplot", "ggpubr",
  "rlang"
)

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

if(length(new_packages) > 0){
  install.packages(new_packages, dep=TRUE)
}

# Load packages
for(package.i in list_of_packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}
rm(list = c("list_of_packages", "new_packages", "package.i"))

### Figures for the UV MR analysis
# Load main MR data
load("#.RData")

main_mr_dat <- lapply(mr_list, `[[`, 1)
main_mr_dat[[7]] <- data.table(id.exposure="Malaria_SMA", id.outcome="Neutrophil", outcome="Neutrophil", exposure="Malaria_SMA", 
  method="Inverse variance weighted", nsnp=NA, b=NA, se=NA, pval=NA, lo_ci=NA, up_ci=NA, or=NA, or_lci95=NA, or_uci95=NA)
main_mr_dat <- rbindlist(main_mr_dat)

main_mr_dat <- main_mr_dat[order(exposure, outcome, method)]
main_mr_dat <- main_mr_dat[method != "Simple mode"]
main_mr_dat <- main_mr_dat[, exposure := fcase(
  exposure == "Malaria_CM", "Cerebral malaria",
  exposure == "Malaria_Other", "Other severe malaria",
  exposure == "Malaria_Overall", "Overall severe malaria",
  exposure == "Malaria_SMA", "Severe malaria anaemia",
  exposure == "Neutrophil", "Neutrophil count"
)][, outcome := fcase(
  outcome == "Malaria_CM", "Cerebral malaria",
  outcome == "Malaria_Other", "Other severe malaria",
  outcome == "Malaria_Overall", "Overall severe malaria",
  outcome == "Malaria_SMA", "Severe malaria anaemia",
  outcome == "Neutrophil", "Neutrophil count"
)]

fwrite(main_mr_dat[method!="MR RAPS", c("exposure", "outcome", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95")], "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")

main_mr_dat <- main_mr_dat[, c("exposure", "outcome", "method", "b", "se", "pval")]


main_mr_dat <- main_mr_dat %>% arrange(factor(main_mr_dat$outcome, 
                             levels=c('Overall severe malaria','Cerebral malaria','Severe malaria anaemia','Other severe malaria','Neutrophil count'),
                             ordered = TRUE))
main_mr_dat$method <- factor(main_mr_dat$method, 
                             levels=rev(c('Wald ratio', 'Inverse variance weighted', 'MR Egger', 'Weighted median', 'Weighted mode', 'MR RAPS')))
main_mr_dat$is_exposure <- ifelse(main_mr_dat$exposure == "Neutrophil count", "Exposure", "Outcome")

# Use Nightingale forestplot package
source("00_ggforestplot.r")

## MR main text
# NC to Malaria, main text
nc_mal_plot <- forestplot_mod(
  df = main_mr_dat[exposure == "Neutrophil count" & method == "Inverse variance weighted", ], name = outcome, estimate = b,
  se = se, pvalue = pval, logodds = TRUE,
  colour = method, shape = method,
  xlab = "Odds ratio for severe malaria (95% CI)
  per 1−SD increment in neutrophil count",
  xlim = c(0.8, 1.2), xtickbreaks = seq(0.0, 2.4, 0.2)
) + theme(axis.text.y = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "MR method", shape = "MR method") + 
  guides(colour = guide_legend(title.position = "top"), shape = guide_legend(title.position = "top"))

main_mr_dat2 <- main_mr_dat %>% arrange(factor(main_mr_dat$exposure, 
                             levels=c('Overall severe malaria','Cerebral malaria','Severe malaria anaemia','Other severe malaria','Neutrophil count'),
                             ordered = TRUE))

# Malaria to NC, main text
mal_nc_plot <- forestplot_mod(
  df = main_mr_dat2[exposure != "Neutrophil count" & method == "Inverse variance weighted", ], name = exposure, estimate = b,
  se = se, pvalue = pval, logodds = FALSE,
  colour = method, shape = method,
  xlab = "SD-unit difference (95% CI) in neutrophil count per
          1 odds ratio increase of severe malaria",
  xlim = c(-1, 2.2), xtickbreaks = seq(0.0, 5.5, 1)
) + theme(axis.text.y=element_blank(),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
  labs(colour = "MR method", shape = "MR method") + 
  guides(colour = guide_legend(title.position = "top"), shape = guide_legend(title.position = "top"))

all_mr_plot <- ggarrange(nc_mal_plot, mal_nc_plot, common.legend = TRUE, labels="AUTO")
ggsave(filename = "#.tiff", 
       plot = all_mr_plot, device="tiff", width = 280, height = 210, dpi=300, compression = "lzw", units="mm", bg="white")

## Load single SNP results
mr_single <- lapply(mr_list, `[[`, 4)
mr_single <- mr_single[-7]
mr_single <- rbindlist(mr_single)
mr_single <- mr_single[, exposure := fcase(
  exposure == "Malaria_CM", "Cerebral malaria",
  exposure == "Malaria_Other", "Other severe malaria",
  exposure == "Malaria_Overall", "Overall severe malaria",
  exposure == "Malaria_SMA", "Severe malaria anaemia",
  exposure == "Neutrophil", "Neutrophil count"
)][, outcome := fcase(
  outcome == "Malaria_CM", "Cerebral malaria",
  outcome == "Malaria_Other", "Other severe malaria",
  outcome == "Malaria_Overall", "Overall severe malaria",
  outcome == "Malaria_SMA", "Severe malaria anaemia",
  outcome == "Neutrophil", "Neutrophil count"
)]
mr_single[, exposure_outcome := paste0(exposure, "_", outcome)]

mr_single_plots <- list()
for(x in 1:length(unique(mr_single$exposure_outcome))) {
  p2 <- mr_forest_plot(mr_single[exposure_outcome == unique(mr_single$exposure_outcome)[x],])
  p2 <- p2[[1]] + labs_pubr() + theme(axis.text.y = element_text(size = 8)) + theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white")
  ) 
  mr_single_plots[[x]] <- p2
 }

mr_nc_plots_combined <- ggarrange(plotlist = mr_single_plots[1:4], labels="AUTO")
ggsave(filename = paste0("#.tiff"), 
  plot = mr_nc_plots_combined, device="tiff", width = 280, height = 210, dpi=300, compression = "lzw", units="mm", bg="white")

mr_mal_plots_combined <- ggarrange(plotlist = mr_single_plots[5:7], labels="AUTO")
ggsave(filename = paste0("#.tiff"), 
  plot = mr_mal_plots_combined, device="tiff", width = 280, height = 210, dpi=300, compression = "lzw", units="mm", bg="white")


## Generate tables with SNPs used in MR analysis
# Neutrophil count table
bolt <- fread("#.bgen")
nc_snps <- mr_single[exposure == "Neutrophil count" & !SNP %in% c("All - Inverse variance weighted", "All - MR Egger"),]
nc_snps <- merge(nc_snps, bolt, by="SNP")
nc_snps <- nc_snps[, c("exposure", "outcome", "SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "b", "se", "p", "BETA", "SE", "P_BOLT_LMM_INF")]
names(nc_snps) <- c("Exposure", "Outcome", "SNP", "CHR", "BP", "NEA", "EA", "EAF", "b.MR", "se.MR", "p.MR", "BETA.BOLT", "SE.BOLT", "P.BOLT")
all_sumstats <- fread("#.txt")
all_sumstats <- all_sumstats[r0.001_clump_lead == "Yes", 1:15][order(CHR, `BP (GRCh37)`)]

nc_snps[, MR_proxy_for := fcase(
  SNP == "rs12758755", "rs12747038", 
  SNP == "rs2325919", "rs2814778",
  SNP == "rs144109344", "rs144109344",
  SNP == "rs73364428", "rs73364428",
  SNP == "rs7460611", "rs10096834", 
  default = "No proxy"
  
)]

fwrite(nc_snps, "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")

# Make tables with cojo index, secondary and MR SNPs
nc_snps_vec <- unique(nc_snps$SNP)
cojo_snps <- fread("#.cojo")
nc_snps_vec <- c(nc_snps_vec, cojo_snps$SNP)
bolt_fuma_snps <- bolt[SNP %in% nc_snps_vec,]

fwrite(bolt_fuma_snps, "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")
fwrite(bolt_fuma_snps[, c("SNP", "CHR", "BP")], "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")

all_sumstats <- fread("#.txt")
plink_bim <- fread("#.bim")
all_sumstats <- all_sumstats[SNP %in% plink_bim$V2, ]
fwrite(all_sumstats, "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")

# Severe malaria table
overall_sm <- fread("#.txt")
cm_sm <- fread("#.txt")
sma_sm <- fread("#.txt")
other_sm <- fread("#.txt")

overall_sm <- overall_sm[!duplicated(MarkerName), ]
cm_sm <- cm_sm[!duplicated(MarkerName), ]
sma_sm <- sma_sm[!duplicated(MarkerName), ]
other_sm <- other_sm[!duplicated(MarkerName), ]

sm_snps <- mr_single[exposure != "Neutrophil count" & !SNP %in% c("All - Inverse variance weighted", "All - MR Egger"),]

overall_sm_snps <- sm_snps[exposure == "Overall severe malaria",]
cm_sm_snps <- sm_snps[exposure == "Cerebral malaria",]
other_sm_snps <- sm_snps[exposure == "Other severe malaria",]


overall_sm_snps <- merge(overall_sm_snps, overall_sm, by.x="SNP", by.y="MarkerName")
cm_sm_snps <- merge(cm_sm_snps, cm_sm, by.x="SNP", by.y="MarkerName")
other_sm_snps <- merge(other_sm_snps, other_sm, by.x="SNP", by.y="MarkerName")

overall_sm_snps <- overall_sm_snps[, c("exposure", "outcome", "SNP", "Chromosome", "Position", "Allele1", "Allele2", "Freq1", "b", "se", "p", "Effect", "StdErr", "P-value")]
names(overall_sm_snps) <- c("Exposure", "Outcome", "SNP", "CHR", "BP", "EA", "NEA", "EAF", "b.MR", "se.MR", "p.MR", "BETA.SNP", "SE.SNP", "P.SNP")

cm_sm_snps <- cm_sm_snps[, c("exposure", "outcome", "SNP", "Chromosome", "Position", "Allele1", "Allele2", "Freq1", "b", "se", "p", "Effect", "StdErr", "P-value")]
names(cm_sm_snps) <- c("Exposure", "Outcome", "SNP", "CHR", "BP", "EA", "NEA", "EAF", "b.MR", "se.MR", "p.MR", "BETA.SNP", "SE.SNP", "P.SNP")

other_sm_snps <- other_sm_snps[, c("exposure", "outcome", "SNP", "Chromosome", "Position", "Allele1", "Allele2", "Freq1", "b", "se", "p", "Effect", "StdErr", "P-value")]
names(other_sm_snps) <- c("Exposure", "Outcome", "SNP", "CHR", "BP", "EA", "NEA", "EAF", "b.MR", "se.MR", "p.MR", "BETA.SNP", "SE.SNP", "P.SNP")

sm_snps <- rbind(overall_sm_snps, cm_sm_snps, other_sm_snps)
sm_snps[, EA := toupper(EA)][, NEA := toupper(NEA)][order(CHR, BP)]

fwrite(sm_snps, "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")

# Make tables with MR SNPs
sm_snps_fuma <- sm_snps[!duplicated(SNP), ]

fwrite(sm_snps_fuma, "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")
fwrite(sm_snps_fuma[, c("SNP", "CHR", "BP")], "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")




