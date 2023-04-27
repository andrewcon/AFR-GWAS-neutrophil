#!/usr/bin/env Rscript

# Load relevant libraries------
list_of_packages <- c(
  "data.table", "tidyverse",
  "ggplot2", "ggpubr", "ggrepel", "RColorBrewer",
  "ggforestplot"
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

# Process GWAS------
# Load all_gwas.cojo.jma file and make a vector of all unique SNPs
all_cojo <- fread("#.cojo")
all_cojo$Chr <- as.character(all_cojo$Chr)
unique_snps <- unique(all_cojo$SNP)
# Create variable with unique_snps to be read with grep in fread
unique_snps_expression <- paste(unique_snps, collapse = "|")

# Pull these out from bolt, all Kpops, metas and Chen paper
# Process chen paper first
chen_gwas <- fread("#.txt")
chen_gwas <- chen_gwas[ID %in% unique_snps,]
chen_gwas$GWAS <- "Chen"
chen_gwas <- chen_gwas[, c("ID", "CHR", "BP", "reference_allele", "other_allele", "eaf", "beta", "se", "p-value", "n_samples", "GWAS")]
setnames(chen_gwas, c("ID", "CHR", "BP", "reference_allele", "other_allele", "eaf", "beta", "se", "p-value", "n_samples", "GWAS"), 
  c("SNP", "CHR", "BP", "EA", "NEA", "eaf", "beta", "se", "pval", "n", "GWAS"))
chen_gwas[SNP == "rs2814778", pval := 1.67E-107]

# Process Kpops
res_wod_K1 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wod_K2 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wod_K3 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wod_K4 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wod_K5 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wod_K6 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wod_K7 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wd_K1 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wd_K2 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wd_K3 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wd_K4 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wd_K5 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wd_K6 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wd_K7 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.txt"), select = c(2:6, 19, 24, 25, 22, 18), header=FALSE)
res_wod_K1$GWAS <- "res_wod_K1"
res_wod_K2$GWAS <- "res_wod_K2"
res_wod_K3$GWAS <- "res_wod_K3"
res_wod_K4$GWAS <- "res_wod_K4"
res_wod_K5$GWAS <- "res_wod_K5"
res_wod_K6$GWAS <- "res_wod_K6"
res_wod_K7$GWAS <- "res_wod_K7"
res_wd_K1$GWAS <- "res_wd_K1"
res_wd_K2$GWAS <- "res_wd_K2"
res_wd_K3$GWAS <- "res_wd_K3"
res_wd_K4$GWAS <- "res_wd_K4"
res_wd_K5$GWAS <- "res_wd_K5"
res_wd_K6$GWAS <- "res_wd_K6"
res_wd_K7$GWAS <- "res_wd_K7"

kpops_all <- rbind(res_wod_K1, res_wod_K2, res_wod_K3, res_wod_K4, res_wod_K5, res_wod_K6, res_wod_K7,
  res_wd_K1, res_wd_K2, res_wd_K3, res_wd_K4, res_wd_K5, res_wd_K6, res_wd_K7)
colnames(kpops_all) <- c("SNP", "CHR", "BP", "EA", "NEA", "eaf", "beta", "se", "pval", "n", "GWAS")
kpops_all$beta <- -1 * kpops_all$beta
kpops_all <- kpops_all[, c("GWAS", "SNP", "CHR", "BP", "EA", "NEA", "eaf", "beta", "se", "pval", "n")]

# Load bolt
bolt <- fread("#.ma")
bolt <- bolt[SNP %in% unique_snps,]
bolt$GWAS <- "bolt"
bolt[, CHR:=NA][, BP:=NA]
bolt <- bolt[, c("GWAS", "SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "P_BOLT_LMM_INF", "N")]
colnames(bolt) <- c("GWAS", "SNP", "CHR", "BP", "EA", "NEA", "eaf", "beta", "se", "pval", "n")

# Load metas
res_wod_meta <- fread("#.ma")
res_wod_meta <- res_wod_meta[rsid %in% unique_snps,]
res_wod_meta$GWAS <- "res_wod_meta"
res_wod_meta[, CHR:=NA][, BP:=NA]
res_wod_meta <- res_wod_meta[, c("GWAS", "rsid", "CHR", "BP", "allele_B", "allele_A", "alleleB_freq", "beta", "se", "P_value", "N")]
colnames(res_wod_meta) <- c("GWAS", "SNP", "CHR", "BP", "EA", "NEA", "eaf", "beta", "se", "pval", "n")

res_wd_meta <- fread("#.ma")
res_wd_meta <- res_wd_meta[rsid %in% unique_snps,]
res_wd_meta$GWAS <- "res_wd_meta"
res_wd_meta[, CHR:=NA][, BP:=NA]
res_wd_meta <- res_wd_meta[, c("GWAS", "rsid", "CHR", "BP", "allele_B", "allele_A", "alleleB_freq", "beta", "se", "P_value", "N")]
colnames(res_wd_meta) <- c("GWAS", "SNP", "CHR", "BP", "EA", "NEA", "eaf", "beta", "se", "pval", "n")

save(all_cojo, chen_gwas, kpops_all, bolt, res_wod_meta, res_wd_meta, file = "#.RData")


# Create a forest plot for all SNP estimates in the relevant datasets
tts_combined <- rbind(chen_gwas, kpops_all, bolt, res_wod_meta, res_wd_meta)
tts_combined[, LowerLimit := beta - se*1.96]
tts_combined[, UpperLimit := beta + se*1.96]
tts_combined[, sig := ifelse(pval < 5e-8, "Yes", "No")]
tts_combined[, GWAS := fcase(
  GWAS == "bolt", "BOLT-LMM",
  GWAS == "Chen", "Chen et al.",
  GWAS == "res_wd_K1", "WD-K1",
  GWAS == "res_wd_K2", "WD-K2",
  GWAS == "res_wd_K3", "WD-K3",
  GWAS == "res_wd_K4", "WD-K4",
  GWAS == "res_wd_K5", "WD-K5",
  GWAS == "res_wd_K6", "WD-K6",
  GWAS == "res_wd_K7", "WD-K7",
  GWAS == "res_wd_meta", "WD-META",
  GWAS == "res_wod_K1", "WOD-K1",
  GWAS == "res_wod_K2", "WOD-K2",
  GWAS == "res_wod_K3", "WOD-K3",
  GWAS == "res_wod_K4", "WOD-K4",
  GWAS == "res_wod_K5", "WOD-K5",
  GWAS == "res_wod_K6", "WOD-K6",
  GWAS == "res_wod_K7", "WOD-K7",
  GWAS == "res_wod_meta", "WOD-META"
)]

all_cojo_shared <- dcast(all_cojo, SNP ~ GWAS, value.var = c("GWAS"))
all_cojo_shared$shared <- paste(all_cojo_shared$bolt, all_cojo_shared$snptest_wd_meta, all_cojo_shared$snptest_wod_meta, sep = "|")
all_cojo_shared <- all_cojo_shared[, c("SNP", "shared")]
all_cojo_shared <- merge(all_cojo_shared, bolt[, c("SNP", "beta", "se", "pval", "eaf", "n")], by = "SNP", all=TRUE)
all_cojo_shared$eac.bolt <- all_cojo_shared$eaf * all_cojo_shared$n
all_cojo_shared <- merge(all_cojo_shared, res_wod_meta[, c("SNP", "beta", "se", "pval", "eaf", "n")], by = "SNP", suffixes = c(".bolt", ".meta_wod"), all=TRUE)
wod_meta_i2 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.unfiltered"), select = c(2, 12))
colnames(wod_meta_i2) <- c("SNP", "i2.meta_wod")
all_cojo_shared <- merge(all_cojo_shared, wod_meta_i2, by="SNP", all=TRUE)
all_cojo_shared$eac.meta_wod <- all_cojo_shared$eaf.meta_wod * all_cojo_shared$n.meta_wod
all_cojo_shared <- merge(all_cojo_shared, res_wd_meta[, c("SNP", "beta", "se", "pval", "eaf", "n")], by = "SNP", all=TRUE)
wd_meta_i2 <- fread(cmd=paste0("grep -wE ", "'",unique_snps_expression, "'", " #.unfiltered"), select = c(2, 12))
colnames(wd_meta_i2) <- c("SNP", "i2.meta_wd")
all_cojo_shared <- merge(all_cojo_shared, wd_meta_i2, by="SNP", all=TRUE)
all_cojo_shared$eac.meta_wd <- all_cojo_shared$eaf * all_cojo_shared$n
colnames(all_cojo_shared)[16:20] <- c("beta.meta_wd", "se.meta_wd", "pval.meta_wd", "eaf.meta_wd", "n.meta_wd")

fwrite(all_cojo_shared, "#.txt", quote=F, row.names=F, sep="\t")

col_wd = c(RColorBrewer::brewer.pal(9, "Blues")[9:3], "grey30")
col_wod = c(RColorBrewer::brewer.pal(9, "Reds")[9:3], "grey60")
pcol = c("black","green4", col_wd, col_wod)
tts_combined <- tts_combined[SNP %in% all_cojo[GWAS == "bolt"]$SNP]

# Create the forest plot
p <- ggplot(data=tts_combined,
    aes(x = GWAS,y = beta, ymin = LowerLimit, ymax = UpperLimit))+
    geom_pointrange(aes(col=GWAS, shape=sig)) + theme_bw() +
    geom_hline(aes(fill=GWAS),yintercept =0, linetype=2)+
    xlab('SNP rsid')+ ylab("Beta (95% Confidence Interval)")+
    geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=GWAS),width=0.25,cex=1) + 
    facet_wrap(~SNP,strip.position="left",nrow=5,scales = "free_y") + 
    theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+ scale_shape_manual(values=c(1, 19) ) +
    coord_flip() + scale_color_manual(values=pcol) +labs(shape="Significant (P<5e-8)")
ggsave(plot = p, 
       filename = "#.tiff", 
       device="tiff", width = 210, height = 280, dpi=300, compression = "lzw", units="mm")

# Create a scatter plot for all SNP estimates between bolt~meta_wd, bolt~meta_wod and meta_wod~meta_wd
p_scatter_1 <- ggscatter(all_cojo_shared, x = "beta.bolt", y = c("beta.meta_wod"), color = "shared",
  combine = FALSE, palette = brewer.pal(5, "Set1"),
  add = "reg.line", xlab="BETA.BOLT", ylab="BETA.META-WOD",
   add.params = list(color = "blue", fill = "lightgray"),
   conf.int = TRUE, 
   cor.coef = TRUE, xlim = c(-0.8, 0.8), xticks.by = 0.3,
   cor.coeff.args = list(method = "spearman", label.sep = "\n") ) + 
  scale_color_manual(labels=c("bolt|NA|NA" = "BOLT", "bolt|NA|snptest_wod_meta"="BOLT & META-WOD", "NA|snptest_wd_meta|NA"="META-WD", "NA|snptest_wd_meta|snptest_wod_meta"="META-WOD & META-WOD"), values = c("red", "blue", "green", "purple"))
p_scatter_2 <- ggscatter(all_cojo_shared, x = "beta.bolt", y = c("beta.meta_wd"), color = "shared",
  combine = FALSE, palette = brewer.pal(5, "Set1"),
  add = "reg.line", xlab="BETA.BOLT", ylab="BETA.META-WD",
   add.params = list(color = "blue", fill = "lightgray"),
   conf.int = TRUE, 
   cor.coef = TRUE, xlim = c(-0.8, 0.8), xticks.by = 0.3,
   cor.coeff.args = list(method = "spearman", label.sep = "\n") ) + 
  scale_color_manual(labels=c("bolt|NA|NA" = "BOLT", "bolt|NA|snptest_wod_meta"="BOLT & META-WOD", "NA|snptest_wd_meta|NA"="META-WD", "NA|snptest_wd_meta|snptest_wod_meta"="META-WOD & META-WOD"), values = c("red", "blue", "green", "purple"))
p_scatter_3 <- ggscatter(all_cojo_shared, x = "beta.meta_wod", y = c("beta.meta_wd"), color = "shared",
  combine = FALSE, palette = brewer.pal(5, "Set1"),
  add = "reg.line", xlab="BETA.META-WOD", ylab="BETA.META-WD",
   add.params = list(color = "blue", fill = "lightgray"),
   conf.int = TRUE, 
   cor.coef = TRUE, xlim = c(-0.8, 0.8), xticks.by = 0.3,
   cor.coeff.args = list(method = "spearman", label.sep = "\n") ) + 
  scale_color_manual(labels=c("bolt|NA|NA" = "BOLT", "bolt|NA|snptest_wod_meta"="BOLT & META-WOD", "NA|snptest_wd_meta|NA"="META-WD", "NA|snptest_wd_meta|snptest_wod_meta"="META-WOD & META-WOD"), values = c("red", "blue", "green", "purple"))
p_scatter_arrange <- ggarrange(p_scatter_1, p_scatter_2, p_scatter_3, ncol = 3, nrow = 1, common.legend = TRUE, align = c("hv"), labels="AUTO")
ggsave(plot = p_scatter_arrange, 
       filename = "#.tiff", 
       device="tiff", width = 240, height = 110, dpi=300, compression = "lzw", units="mm", bg="white")

