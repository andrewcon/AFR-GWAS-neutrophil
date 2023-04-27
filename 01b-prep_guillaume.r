#!/usr/bin/env Rscript

# Load relevant libraries
library(data.table)
library(ieugwasr)
library(glue)

# Load HRC file------
hrc <- fread("#/HRC_GRCh37.vcf", 
             stringsAsFactors = F, skip = 49, select = c("#CHROM", "POS", "ID", "REF", "ALT"))
setnames(hrc, c("#CHROM", "POS", "ID", "REF", "ALT"), c("CHR", "BP", "ID", "HRC_REF", "HRC_ALT"))

# Create vector which will have all SNPs from all immune cells
snp_vec <- c()

in_file <- paste0("#.txt")
mydat <- fread(in_file, stringsAsFactors=F, 
              select = c("rs_number", "reference_allele", "other_allele", 
                        "eaf", "beta", "se", "p-value", "q_statistic", 
                        "q_p-value", "i2", "n_studies", "n_samples", "effects"))
mydat <- mydat[(mydat$'i2' > 0.40 & mydat$'q_p-value' > 0.05) | (mydat$'i2' <= 0.40) | is.nan(mydat$'i2'),]
mydat[, c("CHR", "BP") := tstrsplit(rs_number, ":", fixed=TRUE)]
mydat[, c("BP", "minor_a", "major_a") := tstrsplit(BP, "_", fixed=TRUE)]
mydat$BP <- as.integer(mydat$BP)
mydat <- merge(mydat, hrc, by = c("CHR", "BP"))
mydat[, ID := ifelse(ID == ".", rs_number, ID)]
mydat_exp <- mydat[mydat$'p-value' < 5e-8,]

# Use ieugwasr to clump
setnames(mydat_exp, c("ID", "p-value"), c("rsid", "pval"))
mydat_clump <- ld_clump(dat=mydat_exp, clump_kb=500, clump_r2=0.5, clump_p=0.99, pop = "AFR")
fwrite(mydat_exp, "#.txt", quote=F, row.names=F, sep="\t")
clump_rsids <- c(mydat_clump$rsid, "rs2814778")
writeLines(clump_rsids, "#.txt")

# Get dosage data using PLINK2
system("module load apps/plink/2.0.0")
dos_list <- list()
for (x in 1:7) {
    bgen_file <- paste0("#", x, ".bgen")
    sample_file <- paste0("#", x, ".sample")
    rsid_file <- "#.txt"
    tmp <- tempfile()
    system(glue::glue("plink2 --bgen {bgen_file} --sample {sample_file}",
                    " --extract {rsid_file}",
                    " --export A 'include-alt'",
                    " --threads 8",
                    " --out {tmp}"))
    dos <- data.table::fread(paste0(tmp, ".raw"), data.table = TRUE)
    dos_list[[x]] <- dos
}

# Cbind all chr SNPs into one dosage file
dos_all <- do.call(cbind, dos_list)
dos_all <- dos_all[, -c("IID")]
setnames(dos_all, "FID", "genid")
dos_all <- dos_all[, -c("FID", "PAT", "MAT", "SEX", "PHENOTYPE")]

# Load pheno_afr_final and cbind with it to for further analyses
covars <- fread("#.txt")
covars_dos <- merge(covars, dos_all, by="genid")

# Write this to pheno-cov folder
fwrite(covars_dos, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
