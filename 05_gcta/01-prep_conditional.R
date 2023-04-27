#!/usr/bin/env Rscript

library(data.table)

gwas_paths <- c(
    "#.bgen",
    "#",
    "#"
)
output_paths <- c(
    "#.ma",
    "#.ma",
    "#.ma"
)

bolt <- fread(gwas_paths[1], stringsAsFactors=F)
bolt <- bolt[, c("SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "P_BOLT_LMM_INF")]
bolt$N <- 5976
fwrite(bolt, output_paths[1], quote = F, row.names = F, na = NA, sep = "\t")

snptest_wd_meta <- fread(gwas_paths[2], stringsAsFactors=F)
snptest_wd_meta$alleleB_freq <- 1 - snptest_wd_meta$coded_af
snptest_wd_meta$beta <- -1 * snptest_wd_meta$beta
snptest_wd_meta$N <- snptest_wd_meta$cohort_1_n + snptest_wd_meta$cohort_2_n + snptest_wd_meta$cohort_3_n +
  snptest_wd_meta$cohort_4_n + snptest_wd_meta$cohort_5_n + snptest_wd_meta$cohort_6_n + snptest_wd_meta$cohort_7_n
snptest_wd_meta <- snptest_wd_meta[, c("rsid", "allele_B", "allele_A", "alleleB_freq", "beta", "se", "P_value", "N")]
snptest_wd_meta[, allele_temp := allele_B][, allele_B := allele_A][, allele_A := allele_temp][, allele_temp := NULL]
fwrite(snptest_wd_meta, output_paths[2], quote = F, row.names = F, na = NA, sep = "\t")

snptest_wod_meta <- fread(gwas_paths[3], stringsAsFactors=F)
snptest_wod_meta$alleleB_freq <- 1 - snptest_wod_meta$coded_af
snptest_wod_meta$beta <- -1 * snptest_wod_meta$beta
snptest_wod_meta$N <- snptest_wod_meta$cohort_1_n + snptest_wod_meta$cohort_2_n + snptest_wod_meta$cohort_3_n +
  snptest_wod_meta$cohort_4_n + snptest_wod_meta$cohort_5_n + snptest_wod_meta$cohort_6_n + snptest_wod_meta$cohort_7_n
snptest_wod_meta <- snptest_wod_meta[, c("rsid", "allele_B", "allele_A", "alleleB_freq", "beta", "se", "P_value", "N")]
snptest_wod_meta[, allele_temp := allele_B][, allele_B := allele_A][, allele_A := allele_temp][, allele_temp := NULL]
fwrite(snptest_wod_meta, output_paths[3], quote = F, row.names = F, na = NA, sep = "\t")
