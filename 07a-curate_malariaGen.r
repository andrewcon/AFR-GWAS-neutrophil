# Load packages------
load_pkg <- rlang::quos(tidyverse, data.table, ggplot2, ggpubr, ggrepel, RColorBrewer, car, speedglm, foreach, doParallel, kableExtra, meta, qqman)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

### Load the overall severe malaria data, by population
mgen_all <- fread("#")
keep_cols <- names(mgen_all)[!grepl(c("Vietnam|PNG"), names(mgen_all))]
mgen_all <- mgen_all[, ..keep_cols]

mgen_gambia <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Gambia:N", "Gambia:B_allele_frequency", "Gambia:all_info", "Gambia:beta_1:SM", "Gambia:se_1", "Gambia:pvalue", "Gambia:trusted")]
mgen_gambia[`Gambia:trusted` == 0, `Gambia:beta_1:SM` := NA][`Gambia:trusted` == 0, `Gambia:se_1` := NA][`Gambia:trusted` == 0, `Gambia:pvalue` := NA][,`Gambia:trusted` := NULL]
setnames(mgen_gambia, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_mali <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Mali:N", "Mali:B_allele_frequency", "Mali:all_info", "Mali:beta_1:SM", "Mali:se_1", "Mali:pvalue", "Mali:trusted")]
mgen_mali[`Mali:trusted` == 0, `Mali:beta_1:SM` := NA][`Mali:trusted` == 0, `Mali:se_1` := NA][`Mali:trusted` == 0, `Mali:pvalue` := NA][,`Mali:trusted` := NULL]
setnames(mgen_mali, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_bf <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "BurkinaFaso:N", "BurkinaFaso:B_allele_frequency", "BurkinaFaso:all_info", "BurkinaFaso:beta_1:SM", "BurkinaFaso:se_1", "BurkinaFaso:pvalue", "BurkinaFaso:trusted")]
mgen_bf[`BurkinaFaso:trusted` == 0, `BurkinaFaso:beta_1:SM` := NA][`BurkinaFaso:trusted` == 0, `BurkinaFaso:se_1` := NA][`BurkinaFaso:trusted` == 0, `BurkinaFaso:pvalue` := NA][,`BurkinaFaso:trusted` := NULL]
setnames(mgen_bf, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_ghana <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Ghana:N", "Ghana:B_allele_frequency", "Ghana:all_info", "Ghana:beta_1:SM", "Ghana:se_1", "Ghana:pvalue", "Ghana:trusted")]
mgen_ghana[`Ghana:trusted` == 0, `Ghana:beta_1:SM` := NA][`Ghana:trusted` == 0, `Ghana:se_1` := NA][`Ghana:trusted` == 0, `Ghana:pvalue` := NA][,`Ghana:trusted` := NULL]
setnames(mgen_ghana, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_nigeria <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Nigeria:N", "Nigeria:B_allele_frequency", "Nigeria:all_info", "Nigeria:beta_1:SM", "Nigeria:se_1", "Nigeria:pvalue", "Nigeria:trusted")]
mgen_nigeria[`Nigeria:trusted` == 0, `Nigeria:beta_1:SM` := NA][`Nigeria:trusted` == 0, `Nigeria:se_1` := NA][`Nigeria:trusted` == 0, `Nigeria:pvalue` := NA][,`Nigeria:trusted` := NULL]
setnames(mgen_nigeria, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_cameroon <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Cameroon:N", "Cameroon:B_allele_frequency", "Cameroon:all_info", "Cameroon:beta_1:SM", "Cameroon:se_1", "Cameroon:pvalue", "Cameroon:trusted")]
mgen_cameroon[`Cameroon:trusted` == 0, `Cameroon:beta_1:SM` := NA][`Cameroon:trusted` == 0, `Cameroon:se_1` := NA][`Cameroon:trusted` == 0, `Cameroon:pvalue` := NA][,`Cameroon:trusted` := NULL]
setnames(mgen_cameroon, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_malawi <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Malawi:N", "Malawi:B_allele_frequency", "Malawi:all_info", "Malawi:beta_1:SM", "Malawi:se_1", "Malawi:pvalue", "Malawi:trusted")]
mgen_malawi[`Malawi:trusted` == 0, `Malawi:beta_1:SM` := NA][`Malawi:trusted` == 0, `Malawi:se_1` := NA][`Malawi:trusted` == 0, `Malawi:pvalue` := NA][,`Malawi:trusted` := NULL]
setnames(mgen_malawi, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_tanzania <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Tanzania:N", "Tanzania:B_allele_frequency", "Tanzania:all_info", "Tanzania:beta_1:SM", "Tanzania:se_1", "Tanzania:pvalue", "Tanzania:trusted")]
mgen_tanzania[`Tanzania:trusted` == 0, `Tanzania:beta_1:SM` := NA][`Tanzania:trusted` == 0, `Tanzania:se_1` := NA][`Tanzania:trusted` == 0, `Tanzania:pvalue` := NA][,`Tanzania:trusted` := NULL]
setnames(mgen_tanzania, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_kenya <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Kenya:N", "Kenya:B_allele_frequency", "Kenya:all_info", "Kenya:beta_1:SM", "Kenya:se_1", "Kenya:pvalue", "Kenya:trusted")]
mgen_kenya[`Kenya:trusted` == 0, `Kenya:beta_1:SM` := NA][`Kenya:trusted` == 0, `Kenya:se_1` := NA][`Kenya:trusted` == 0, `Kenya:pvalue` := NA][,`Kenya:trusted` := NULL]
setnames(mgen_kenya, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

rm(mgen_all)

# Write to file
fwrite(mgen_gambia, "#/Gambia.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_mali, "#/Mali.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_bf, "#/BurkinaFaso.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_ghana, "#/Ghana.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_nigeria, "#/Nigeria.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_cameroon, "#/Cameroon.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_malawi, "#/Malawi.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_tanzania, "#/Tanzania.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_kenya, "#/Kenya.txt", row.names=F, na=NA, quote=F, sep="\t")


### Load the sub-phenotype files, by population
mgen_all <- fread("#")
keep_cols <- names(mgen_all)[!grepl(c("Vietnam|PNG"), names(mgen_all))]
mgen_all <- mgen_all[, ..keep_cols]

## Curating CM
mgen_gambia <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Gambia:N", "Gambia:B_allele_frequency", "Gambia:all_info", "Gambia:beta_1:CM", "Gambia:se_1", "Gambia:pvalue", "Gambia:trusted")]
mgen_gambia[`Gambia:trusted` == 0, `Gambia:beta_1:CM` := NA][`Gambia:trusted` == 0, `Gambia:se_1` := NA][`Gambia:trusted` == 0, `Gambia:pvalue` := NA][,`Gambia:trusted` := NULL]
setnames(mgen_gambia, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_mali <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Mali:N", "Mali:B_allele_frequency", "Mali:all_info", "Mali:beta_1:CM", "Mali:se_1", "Mali:pvalue", "Mali:trusted")]
mgen_mali[`Mali:trusted` == 0, `Mali:beta_1:CM` := NA][`Mali:trusted` == 0, `Mali:se_1` := NA][`Mali:trusted` == 0, `Mali:pvalue` := NA][,`Mali:trusted` := NULL]
setnames(mgen_mali, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_bf <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "BurkinaFaso:N", "BurkinaFaso:B_allele_frequency", "BurkinaFaso:all_info", "BurkinaFaso:beta_1:CM", "BurkinaFaso:se_1", "BurkinaFaso:pvalue", "BurkinaFaso:trusted")]
mgen_bf[`BurkinaFaso:trusted` == 0, `BurkinaFaso:beta_1:CM` := NA][`BurkinaFaso:trusted` == 0, `BurkinaFaso:se_1` := NA][`BurkinaFaso:trusted` == 0, `BurkinaFaso:pvalue` := NA][,`BurkinaFaso:trusted` := NULL]
setnames(mgen_bf, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_ghana <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Ghana:N", "Ghana:B_allele_frequency", "Ghana:all_info", "Ghana:beta_1:CM", "Ghana:se_1", "Ghana:pvalue", "Ghana:trusted")]
mgen_ghana[`Ghana:trusted` == 0, `Ghana:beta_1:CM` := NA][`Ghana:trusted` == 0, `Ghana:se_1` := NA][`Ghana:trusted` == 0, `Ghana:pvalue` := NA][,`Ghana:trusted` := NULL]
setnames(mgen_ghana, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_nigeria <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Nigeria:N", "Nigeria:B_allele_frequency", "Nigeria:all_info", "Nigeria:beta_1:CM", "Nigeria:se_1", "Nigeria:pvalue", "Nigeria:trusted")]
mgen_nigeria[`Nigeria:trusted` == 0, `Nigeria:beta_1:CM` := NA][`Nigeria:trusted` == 0, `Nigeria:se_1` := NA][`Nigeria:trusted` == 0, `Nigeria:pvalue` := NA][,`Nigeria:trusted` := NULL]
setnames(mgen_nigeria, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_cameroon <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Cameroon:N", "Cameroon:B_allele_frequency", "Cameroon:all_info", "Cameroon:beta_1:CM", "Cameroon:se_1", "Cameroon:pvalue", "Cameroon:trusted")]
mgen_cameroon[`Cameroon:trusted` == 0, `Cameroon:beta_1:CM` := NA][`Cameroon:trusted` == 0, `Cameroon:se_1` := NA][`Cameroon:trusted` == 0, `Cameroon:pvalue` := NA][,`Cameroon:trusted` := NULL]
setnames(mgen_cameroon, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_malawi <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Malawi:N", "Malawi:B_allele_frequency", "Malawi:all_info", "Malawi:beta_1:CM", "Malawi:se_1", "Malawi:pvalue", "Malawi:trusted")]
mgen_malawi[`Malawi:trusted` == 0, `Malawi:beta_1:CM` := NA][`Malawi:trusted` == 0, `Malawi:se_1` := NA][`Malawi:trusted` == 0, `Malawi:pvalue` := NA][,`Malawi:trusted` := NULL]
setnames(mgen_malawi, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_tanzania <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Tanzania:N", "Tanzania:B_allele_frequency", "Tanzania:all_info", "Tanzania:beta_1:CM", "Tanzania:se_1", "Tanzania:pvalue", "Tanzania:trusted")]
mgen_tanzania[`Tanzania:trusted` == 0, `Tanzania:beta_1:CM` := NA][`Tanzania:trusted` == 0, `Tanzania:se_1` := NA][`Tanzania:trusted` == 0, `Tanzania:pvalue` := NA][,`Tanzania:trusted` := NULL]
setnames(mgen_tanzania, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_kenya <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Kenya:N", "Kenya:B_allele_frequency", "Kenya:all_info", "Kenya:beta_1:CM", "Kenya:se_1", "Kenya:pvalue", "Kenya:trusted")]
mgen_kenya[`Kenya:trusted` == 0, `Kenya:beta_1:CM` := NA][`Kenya:trusted` == 0, `Kenya:se_1` := NA][`Kenya:trusted` == 0, `Kenya:pvalue` := NA][,`Kenya:trusted` := NULL]
setnames(mgen_kenya, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

# Write to file
fwrite(mgen_gambia, "#/CM_Gambia.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_mali, "#/CM_Mali.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_bf, "#/CM_BurkinaFaso.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_ghana, "#/CM_Ghana.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_nigeria, "#/CM_Nigeria.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_cameroon, "#/CM_Cameroon.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_malawi, "#/CM_Malawi.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_tanzania, "#/CM_Tanzania.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_kenya, "#/CM_Kenya.txt", row.names=F, na=NA, quote=F, sep="\t")


## Curating SMA
mgen_gambia <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Gambia:N", "Gambia:B_allele_frequency", "Gambia:all_info", "Gambia:beta_3:SMA", "Gambia:se_3", "Gambia:pvalue", "Gambia:trusted")]
mgen_gambia[`Gambia:trusted` == 0, `Gambia:beta_3:SMA` := NA][`Gambia:trusted` == 0, `Gambia:se_3` := NA][`Gambia:trusted` == 0, `Gambia:pvalue` := NA][,`Gambia:trusted` := NULL]
setnames(mgen_gambia, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_mali <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Mali:N", "Mali:B_allele_frequency", "Mali:all_info", "Mali:beta_3:SMA", "Mali:se_3", "Mali:pvalue", "Mali:trusted")]
mgen_mali[`Mali:trusted` == 0, `Mali:beta_3:SMA` := NA][`Mali:trusted` == 0, `Mali:se_3` := NA][`Mali:trusted` == 0, `Mali:pvalue` := NA][,`Mali:trusted` := NULL]
setnames(mgen_mali, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_bf <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "BurkinaFaso:N", "BurkinaFaso:B_allele_frequency", "BurkinaFaso:all_info", "BurkinaFaso:beta_3:SMA", "BurkinaFaso:se_3", "BurkinaFaso:pvalue", "BurkinaFaso:trusted")]
mgen_bf[`BurkinaFaso:trusted` == 0, `BurkinaFaso:beta_3:SMA` := NA][`BurkinaFaso:trusted` == 0, `BurkinaFaso:se_3` := NA][`BurkinaFaso:trusted` == 0, `BurkinaFaso:pvalue` := NA][,`BurkinaFaso:trusted` := NULL]
setnames(mgen_bf, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_ghana <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Ghana:N", "Ghana:B_allele_frequency", "Ghana:all_info", "Ghana:beta_3:SMA", "Ghana:se_3", "Ghana:pvalue", "Ghana:trusted")]
mgen_ghana[`Ghana:trusted` == 0, `Ghana:beta_3:SMA` := NA][`Ghana:trusted` == 0, `Ghana:se_3` := NA][`Ghana:trusted` == 0, `Ghana:pvalue` := NA][,`Ghana:trusted` := NULL]
setnames(mgen_ghana, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_nigeria <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Nigeria:N", "Nigeria:B_allele_frequency", "Nigeria:all_info", "Nigeria:beta_3:SMA", "Nigeria:se_3", "Nigeria:pvalue", "Nigeria:trusted")]
mgen_nigeria[`Nigeria:trusted` == 0, `Nigeria:beta_3:SMA` := NA][`Nigeria:trusted` == 0, `Nigeria:se_3` := NA][`Nigeria:trusted` == 0, `Nigeria:pvalue` := NA][,`Nigeria:trusted` := NULL]
setnames(mgen_nigeria, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_cameroon <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Cameroon:N", "Cameroon:B_allele_frequency", "Cameroon:all_info", "Cameroon:beta_3:SMA", "Cameroon:se_3", "Cameroon:pvalue", "Cameroon:trusted")]
mgen_cameroon[`Cameroon:trusted` == 0, `Cameroon:beta_3:SMA` := NA][`Cameroon:trusted` == 0, `Cameroon:se_3` := NA][`Cameroon:trusted` == 0, `Cameroon:pvalue` := NA][,`Cameroon:trusted` := NULL]
setnames(mgen_cameroon, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_malawi <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Malawi:N", "Malawi:B_allele_frequency", "Malawi:all_info", "Malawi:beta_3:SMA", "Malawi:se_3", "Malawi:pvalue", "Malawi:trusted")]
mgen_malawi[`Malawi:trusted` == 0, `Malawi:beta_3:SMA` := NA][`Malawi:trusted` == 0, `Malawi:se_3` := NA][`Malawi:trusted` == 0, `Malawi:pvalue` := NA][,`Malawi:trusted` := NULL]
setnames(mgen_malawi, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_tanzania <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Tanzania:N", "Tanzania:B_allele_frequency", "Tanzania:all_info", "Tanzania:beta_3:SMA", "Tanzania:se_3", "Tanzania:pvalue", "Tanzania:trusted")]
mgen_tanzania[`Tanzania:trusted` == 0, `Tanzania:beta_3:SMA` := NA][`Tanzania:trusted` == 0, `Tanzania:se_3` := NA][`Tanzania:trusted` == 0, `Tanzania:pvalue` := NA][,`Tanzania:trusted` := NULL]
setnames(mgen_tanzania, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_kenya <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Kenya:N", "Kenya:B_allele_frequency", "Kenya:all_info", "Kenya:beta_3:SMA", "Kenya:se_3", "Kenya:pvalue", "Kenya:trusted")]
mgen_kenya[`Kenya:trusted` == 0, `Kenya:beta_3:SMA` := NA][`Kenya:trusted` == 0, `Kenya:se_3` := NA][`Kenya:trusted` == 0, `Kenya:pvalue` := NA][,`Kenya:trusted` := NULL]
setnames(mgen_kenya, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

# Write to file
fwrite(mgen_gambia, "#/SMA_Gambia.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_mali, "#/SMA_Mali.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_bf, "#/SMA_BurkinaFaso.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_ghana, "#/SMA_Ghana.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_nigeria, "#/SMA_Nigeria.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_cameroon, "#/SMA_Cameroon.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_malawi, "#/SMA_Malawi.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_tanzania, "#/SMA_Tanzania.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_kenya, "#/SMA_Kenya.txt", row.names=F, na=NA, quote=F, sep="\t")


## Curating OTHER severe malaria
mgen_gambia <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Gambia:N", "Gambia:B_allele_frequency", "Gambia:all_info", "Gambia:beta_2:OTHER", "Gambia:se_2", "Gambia:pvalue", "Gambia:trusted")]
mgen_gambia[`Gambia:trusted` == 0, `Gambia:beta_2:OTHER` := NA][`Gambia:trusted` == 0, `Gambia:se_2` := NA][`Gambia:trusted` == 0, `Gambia:pvalue` := NA][,`Gambia:trusted` := NULL]
setnames(mgen_gambia, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_mali <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Mali:N", "Mali:B_allele_frequency", "Mali:all_info", "Mali:beta_2:OTHER", "Mali:se_2", "Mali:pvalue", "Mali:trusted")]
mgen_mali[`Mali:trusted` == 0, `Mali:beta_2:OTHER` := NA][`Mali:trusted` == 0, `Mali:se_2` := NA][`Mali:trusted` == 0, `Mali:pvalue` := NA][,`Mali:trusted` := NULL]
setnames(mgen_mali, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_bf <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "BurkinaFaso:N", "BurkinaFaso:B_allele_frequency", "BurkinaFaso:all_info", "BurkinaFaso:beta_2:OTHER", "BurkinaFaso:se_2", "BurkinaFaso:pvalue", "BurkinaFaso:trusted")]
mgen_bf[`BurkinaFaso:trusted` == 0, `BurkinaFaso:beta_2:OTHER` := NA][`BurkinaFaso:trusted` == 0, `BurkinaFaso:se_2` := NA][`BurkinaFaso:trusted` == 0, `BurkinaFaso:pvalue` := NA][,`BurkinaFaso:trusted` := NULL]
setnames(mgen_bf, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_ghana <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Ghana:N", "Ghana:B_allele_frequency", "Ghana:all_info", "Ghana:beta_2:OTHER", "Ghana:se_2", "Ghana:pvalue", "Ghana:trusted")]
mgen_ghana[`Ghana:trusted` == 0, `Ghana:beta_2:OTHER` := NA][`Ghana:trusted` == 0, `Ghana:se_2` := NA][`Ghana:trusted` == 0, `Ghana:pvalue` := NA][,`Ghana:trusted` := NULL]
setnames(mgen_ghana, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_nigeria <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Nigeria:N", "Nigeria:B_allele_frequency", "Nigeria:all_info", "Nigeria:beta_2:OTHER", "Nigeria:se_2", "Nigeria:pvalue", "Nigeria:trusted")]
mgen_nigeria[`Nigeria:trusted` == 0, `Nigeria:beta_2:OTHER` := NA][`Nigeria:trusted` == 0, `Nigeria:se_2` := NA][`Nigeria:trusted` == 0, `Nigeria:pvalue` := NA][,`Nigeria:trusted` := NULL]
setnames(mgen_nigeria, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_cameroon <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Cameroon:N", "Cameroon:B_allele_frequency", "Cameroon:all_info", "Cameroon:beta_2:OTHER", "Cameroon:se_2", "Cameroon:pvalue", "Cameroon:trusted")]
mgen_cameroon[`Cameroon:trusted` == 0, `Cameroon:beta_2:OTHER` := NA][`Cameroon:trusted` == 0, `Cameroon:se_2` := NA][`Cameroon:trusted` == 0, `Cameroon:pvalue` := NA][,`Cameroon:trusted` := NULL]
setnames(mgen_cameroon, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_malawi <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Malawi:N", "Malawi:B_allele_frequency", "Malawi:all_info", "Malawi:beta_2:OTHER", "Malawi:se_2", "Malawi:pvalue", "Malawi:trusted")]
mgen_malawi[`Malawi:trusted` == 0, `Malawi:beta_2:OTHER` := NA][`Malawi:trusted` == 0, `Malawi:se_2` := NA][`Malawi:trusted` == 0, `Malawi:pvalue` := NA][,`Malawi:trusted` := NULL]
setnames(mgen_malawi, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_tanzania <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Tanzania:N", "Tanzania:B_allele_frequency", "Tanzania:all_info", "Tanzania:beta_2:OTHER", "Tanzania:se_2", "Tanzania:pvalue", "Tanzania:trusted")]
mgen_tanzania[`Tanzania:trusted` == 0, `Tanzania:beta_2:OTHER` := NA][`Tanzania:trusted` == 0, `Tanzania:se_2` := NA][`Tanzania:trusted` == 0, `Tanzania:pvalue` := NA][,`Tanzania:trusted` := NULL]
setnames(mgen_tanzania, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

mgen_kenya <- mgen_all[, c("rsid", "chromosome", "position", "alleleA", "alleleB", "Kenya:N", "Kenya:B_allele_frequency", "Kenya:all_info", "Kenya:beta_2:OTHER", "Kenya:se_2", "Kenya:pvalue", "Kenya:trusted")]
mgen_kenya[`Kenya:trusted` == 0, `Kenya:beta_2:OTHER` := NA][`Kenya:trusted` == 0, `Kenya:se_2` := NA][`Kenya:trusted` == 0, `Kenya:pvalue` := NA][,`Kenya:trusted` := NULL]
setnames(mgen_kenya, c("rsid", "chromosome", "position", "alleleA", "alleleB", "all_total", "all_maf", "frequentist_add_info", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue"))

# Write to file
fwrite(mgen_gambia, "#/OTHER_Gambia.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_mali, "#/OTHER_Mali.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_bf, "#/OTHER_BurkinaFaso.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_ghana, "#/OTHER_Ghana.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_nigeria, "#/OTHER_Nigeria.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_cameroon, "#/OTHER_Cameroon.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_malawi, "#/OTHER_Malawi.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_tanzania, "#/OTHER_Tanzania.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_kenya, "#/OTHER_Kenya.txt", row.names=F, na=NA, quote=F, sep="\t")

# Write test
fwrite(mgen_gambia[chromosome=="22",], "#/TEST_OTHER_Gambia.txt", row.names=F, na=NA, quote=F, sep="\t")
fwrite(mgen_mali[chromosome=="22",], "#/TEST_OTHER_Mali.txt", row.names=F, na=NA, quote=F, sep="\t")

