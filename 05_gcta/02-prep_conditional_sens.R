#!/usr/bin/env Rscript

library(data.table)

gwas_paths <- c(
    "#.bgen.gz"
)
output_paths <- c(
    "#.ma"
)

bolt <- fread(gwas_paths[1], stringsAsFactors=F)
bolt <- bolt[, c("SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "P_BOLT_LMM_INF")]
bolt$N <- 5310
fwrite(bolt, output_paths[1], quote = F, row.names = F, na = NA, sep = "\t")
