#!/usr/bin/env Rscript

library(data.table)

### 1) Make the pheno, covar, qcovar files for the GWAS-parameters GREML
ukbb <- fread("#.txt", stringsAsFactors = F)
ukbb[, FID := genid][, IID := genid][, genid := NULL]
ukbb <- ukbb[, c(ncol(ukbb)-1, ncol(ukbb), 1:(ncol(ukbb)-2)), with=F]

# Remove NA for neutrophil count
ukbb <- ukbb[!is.na(nc_log)]

pheno_cols <- c("FID", "IID", "nc_log")
covar_cols <- c("FID", "IID", "sex", "neutrophil_device", "sample_year", "sample_month")
qcovar_cols <- c("FID", "IID", "age", "sample_day", "sample_day_minutes", paste0("PC", 1:100))

pheno <- ukbb[, ..pheno_cols]
covar <- ukbb[, ..covar_cols]
qcovar <- ukbb[, ..qcovar_cols]

fwrite(pheno, "#.pheno", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)
fwrite(covar, "#.covar", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)
fwrite(qcovar, "#.qcovar", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)

### 2) Make the pheno, covar, qcovar files for the GWAS-parameters GREML
ukbb <- fread("#.txt", stringsAsFactors = F)
ukbb[, FID := genid][, IID := genid][, genid := NULL]
ukbb <- ukbb[, c(ncol(ukbb)-1, ncol(ukbb), 1:(ncol(ukbb)-2)), with=F]

# Remove NA for neutrophil count
ukbb <- ukbb[!is.na(nc_log)]

pheno_cols <- c("FID", "IID", "nc_log")
covar_cols <- c("FID", "IID", "sex", "neutrophil_device", "sample_year", "sample_month")
qcovar_cols <- c("FID", "IID", "age", "sample_day", "sample_day_minutes", paste0("PC", 1:100), "rs2814778_T")

pheno <- ukbb[, ..pheno_cols]
covar <- ukbb[, ..covar_cols]
qcovar <- ukbb[, ..qcovar_cols]

fwrite(pheno, "#.pheno", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)
fwrite(covar, "#.covar", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)
fwrite(qcovar, "#.qcovar", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)

### 3) Make the pheno, covar, qcovar AND gxe files for also studying other variables
ukbb <- fread("#.txt", stringsAsFactors = F)
ukbb[, FID := genid][, IID := genid][, genid := NULL]
ukbb <- ukbb[, c(ncol(ukbb)-1, ncol(ukbb), 1:(ncol(ukbb)-2)), with=F]

ukbb <- ukbb[, c("FID", "IID", "nc_log", "neutrophil_device", 
                     "sample_year", "sample_month", "sample_day", "sample_day_minutes", 
                     "sex", "age", "Kpop", 
                     paste0("PC", 1:100),
                     "BMI", "height", "weight", "home_north", 
                     "home_east", "townsendDI", "smoking_status", "drinker_status", 
                     "rs2814778_T")]

# Remove NA for ALL columns
ukbb <- na.omit(ukbb)

pheno_cols <- c("FID", "IID", "nc_log")
covar_cols <- c("FID", "IID", "sex", "neutrophil_device", "sample_year", "sample_month", "Kpop")
qcovar_cols <- c("FID", "IID", "age", "sample_month", "sample_day", "sample_day_minutes", paste0("PC", 1:100))
gxe_cols <- c("FID", "IID", "sex", "BMI", "smoking_status", "drinker_status")

pheno <- ukbb[, ..pheno_cols]
covar <- ukbb[, ..covar_cols]
qcovar <- ukbb[, ..qcovar_cols]
gxe <- ukbb[, ..gxe_cols][, BMI := fcase(
    BMI < 18.5, "Underweight",
    BMI >= 18.5 & BMI < 25, "NormalWeight",
    BMI >= 25 & BMI < 30, "Overweight",
    BMI >= 30 & BMI < 40, "Obese",
    BMI >= 40, "MorbidlyObese"
)]

fwrite(pheno, "#.pheno", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)
fwrite(covar, "#.covar", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)
fwrite(qcovar, "#.qcovar", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)
fwrite(gxe, "#.gxe", row.names = F, col.names = F, quote = F, sep = "\t", na = NA)

### Finally, make a file with all GRM path names
grm_files <- tools::file_path_sans_ext(list.files("#", full.names=TRUE, pattern="*.log"))
grm_files <- grm_files[grm_files != "#"]
writeLines(grm_files, "#")
