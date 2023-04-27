#!/usr/bin/env Rscript

# Load relevant libraries
library(data.table)
library(ieugwasr)
library(glue)

# Load BOLT, meta_wod and meta_wd------
bolt <- fread("#.bgen")
bolt <- bolt[P_BOLT_LMM_INF < 5e-8,]
fwrite(bolt, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
bolt <- fread("#.txt")

meta_wd <- fread("#.txt.unfiltered")
meta_wd$beta <- -1 * meta_wd$beta
meta_wd$coded_af <- 1 - meta_wd$coded_af
meta_wod <- fread("#.txt.unfiltered")
meta_wod$beta <- -1 * meta_wod$beta
meta_wod$coded_af <- 1 - meta_wod$coded_af
fwrite(meta_wd, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(meta_wod, "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Clump BOLT, meta_wod and meta_wd to PRS standards
# r2 = 0.6, kb = 250
system("module load apps/plink/1.90")
bed_path <- paste0("#/data.chr1-22")

# bolt
sumstat_file <- "#.txt"
tmp <- tempfile()
system(glue::glue("plink --bfile {bed_path}",
                " --clump {sumstat_file}",
                " --clump-p1 1",
                " --clump-p2 1",
                " --clump-r2 0.6",
                " --clump-kb 250",
                " --clump-snp-field 'SNP'",
                " --clump-field 'P_BOLT_LMM_INF'",
                " --threads 8",
                " --out {tmp}"))
clump <- data.table::fread(paste0(tmp, ".clumped"), data.table = TRUE, fill=TRUE)
bolt[, r0.6_clump_lead := ifelse(SNP %in% clump$SNP, "Yes", "No")]

# Clump BOLT, meta_wod and meta_wd to PRS standards
# r2 = 0.1, kb = 250
system("module load apps/plink/1.90")
bed_path <- paste0("#")

# bolt
sumstat_file <- "#.txt"
tmp <- tempfile()
system(glue::glue("plink --bfile {bed_path}",
                " --clump {sumstat_file}",
                " --clump-p1 1",
                " --clump-p2 1",
                " --clump-r2 0.1",
                " --clump-kb 250",
                " --clump-snp-field 'SNP'",
                " --clump-field 'P_BOLT_LMM_INF'",
                " --threads 8",
                " --out {tmp}"))
clump <- data.table::fread(paste0(tmp, ".clumped"), data.table = TRUE, fill=TRUE)
bolt[, r0.1_clump_lead := ifelse(SNP %in% clump$SNP, "Yes", "No")]


# meta_wd
sumstat_file <- "#.txt"
tmp <- tempfile()
system(glue::glue("plink --bfile {bed_path}",
                " --clump {sumstat_file}",
                " --clump-p1 1",
                " --clump-p2 1",
                " --clump-r2 0.1",
                " --clump-kb 250",
                " --clump-snp-field 'rsid'",
                " --clump-field 'P_value'",
                " --threads 8",
                " --out {tmp}"))
clump <- data.table::fread(paste0(tmp, ".clumped"), data.table = TRUE, fill=TRUE)
meta_wd[, r0.1_clump_lead := ifelse(rsid %in% clump$SNP, "Yes", "No")]

# meta_wd
sumstat_file <- "#.txt"
tmp <- tempfile()
system(glue::glue("plink --bfile {bed_path}",
                " --clump {sumstat_file}",
                " --clump-p1 1",
                " --clump-p2 1",
                " --clump-r2 0.1",
                " --clump-kb 250",
                " --clump-snp-field 'rsid'",
                " --clump-field 'P_value'",
                " --threads 8",
                " --out {tmp}"))
clump <- data.table::fread(paste0(tmp, ".clumped"), data.table = TRUE, fill=TRUE)
meta_wod[, r0.1_clump_lead := ifelse(rsid %in% clump$SNP, "Yes", "No")]

# Clump BOLT, meta_wod and meta_wd to MR standards
# r2 = 0.001, kb = 10000
# bolt
sumstat_file <- "#.txt"
tmp <- tempfile()
system(glue::glue("plink --bfile {bed_path}",
                " --clump {sumstat_file}",
                " --clump-p1 1",
                " --clump-p2 1",
                " --clump-r2 0.001",
                " --clump-kb 10000",
                " --clump-snp-field 'SNP'",
                " --clump-field 'P_BOLT_LMM_INF'",
                " --threads 8",
                " --out {tmp}"))
clump <- data.table::fread(paste0(tmp, ".clumped"), data.table = TRUE, fill=TRUE)
bolt[, r0.001_clump_lead := ifelse(SNP %in% clump$SNP, "Yes", "No")]

# meta_wd
sumstat_file <- "#.txt"
tmp <- tempfile()
system(glue::glue("plink --bfile {bed_path}",
                " --clump {sumstat_file}",
                " --clump-p1 1",
                " --clump-p2 1",
                " --clump-r2 0.001",
                " --clump-kb 10000",
                " --clump-snp-field 'rsid'",
                " --clump-field 'P_value'",
                " --threads 8",
                " --out {tmp}"))
clump <- data.table::fread(paste0(tmp, ".clumped"), data.table = TRUE, fill=TRUE)
meta_wd[, r0.001_clump_lead := ifelse(rsid %in% clump$SNP, "Yes", "No")]

# meta_wd
sumstat_file <- "#.txt"
tmp <- tempfile()
system(glue::glue("plink --bfile {bed_path}",
                " --clump {sumstat_file}",
                " --clump-p1 1",
                " --clump-p2 1",
                " --clump-r2 0.001",
                " --clump-kb 10000",
                " --clump-snp-field 'rsid'",
                " --clump-field 'P_value'",
                " --threads 8",
                " --out {tmp}"))
clump <- data.table::fread(paste0(tmp, ".clumped"), data.table = TRUE, fill=TRUE)
meta_wod[, r0.001_clump_lead := ifelse(rsid %in% clump$SNP, "Yes", "No")]

# Combine with cojo index SNPs
cojo_index <- fread("#.cojo")
cojo_index <- cojo_index[p < 5e-8,]
bolt[, cojo_index := ifelse(SNP %in% cojo_index$SNP, "Yes", "No")]

fwrite(bolt, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(meta_wd, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(meta_wod, "#.txt", quote=F, row.names=F, na=NA, sep="\t")


