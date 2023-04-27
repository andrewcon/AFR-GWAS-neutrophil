# Load packages------
load_pkg <- rlang::quos(ggplot2, ggpubr, ggrepel, RColorBrewer, car, 
    speedglm, foreach, doParallel, kableExtra, meta, qqman, officer, flextable, gtsummary, rlang, 
    openxlsx, readxl, rsnps, biomaRt, splitstackshape, tidyverse, data.table)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

#####
##### 1) Generate table detailing the sample
#####
# Make docx table
FitFlextableToPage <- function(ft, pgwidth = 6){

  ft_out <- ft %>% autofit()

  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}
set_flextable_defaults(
  font.size = 11, font.family = "Arial",
  padding = 4)

# Load covars
covars <- fread("#.txt")
covars$menopause <- as.character(covars$menopause)
covars[is.na(menopause), menopause := "male"]

# Convert to factors
covars$sample_year <- factor(covars$sample_year, 
                             levels = c(2007, 2008, 2009, 2010, 2015), 
                             labels = c(2007, 2008, 2009, 2010, 2015)
)
covars$sample_month <- factor(covars$sample_month, 
                              levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
                              labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
)
covars$sex <- factor(covars$sex, 
                     levels = c(0, 1), 
                     labels = c("Female", "Male")
)
covars$menopause <- factor(covars$menopause, 
                           levels = c("male", "-3", "0", "2", "3", "1"), 
                           labels = c("Male", "Prefer not answer", "No", "Hysterectomy", "Not sure - other", "Yes")
)
covars$cob_un_regions <- as.factor(covars$cob_un_regions)
covars$Kpop <- as.factor(covars$Kpop)
covars$rm_astle <- factor(covars$rm_astle, 
                          levels = c(0, 1), 
                          labels = c("No", "Yes")
)
covars$smoking_status <- factor(covars$smoking_status, 
                                levels = c("prefer_not_answer", "never", "previous", "current"), 
                                labels = c("Prefer not to answer", "Never", "Previous", "Current")
)
covars$drinker_status <- factor(covars$drinker_status, 
                                levels = c("prefer_not_answer", "never", "previous", "current"), 
                                labels = c("Prefer not to answer", "Never", "Previous", "Current")
)
covars$rs2814778_T <- factor(covars$rs2814778_T, 
                             levels = c(0, 1, 2), 
                             labels = c("TT", "TC", "CC")
)

# Generate table and write to .docx
covars_tb <- covars %>% 
  dplyr::select(neutrophil_count, sex, age, BMI, Kpop) %>%
  drop_na(neutrophil_count) %>%
  tbl_summary(     
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = all_categorical() ~ "categorical",                 # force all categorical levels to display
    label  = list(                                              # display labels for column names
      neutrophil_count ~ "Neutrophil count (10^9 cells/Litre)",
      sex   ~ "Genetic sex",                           
      age ~ "Age (years)",
      BMI ~ "BMI (kg/m^2)",
      Kpop    ~ "K-means cluster"),
    missing_text = "Missing"                                    # how missing values should display
  ) %>% as_flex_table() %>% 
    autofit() %>% FitFlextableToPage()

save_as_docx(covars_tb, path = "#.docx")

covars_tb_all <- covars %>% 
  dplyr::select(neutrophil_count, sex, menopause, age, BMI, smoking_status, drinker_status, Kpop) %>%
  drop_na(neutrophil_count) %>%
  tbl_summary(     
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = all_categorical() ~ "categorical",                 # force all categorical levels to display
    label  = list(                                              # display labels for column names
      neutrophil_count ~ "Neutrophil count (10^9 cells/Litre)",
      sex   ~ "Genetic sex",
      menopause ~ "Menopause status",                           
      age ~ "Age (years)",
      BMI ~ "BMI (kg/m^2)",
      smoking_status   ~ "Smoking status",                           
      drinker_status ~ "Alcohol drinker status",
      Kpop    ~ "K-means cluster"),
    missing_text = "Missing"                                    # how missing values should display
  ) %>% as_flex_table() %>% 
    autofit() %>% FitFlextableToPage()

save_as_docx(covars_tb_all, path = "#.docx")


#####
##### 2) make a supplementary table for ALL of the SNPs that Chen and Reiner (PMID: 21738479; 1 more here at CXCL2) 
##### identifed and report all of the GWAS sum stats for each k-pop, meta, and bolt. N, Ef_allele, alt_allele, INFO, EAF, MAC, BETA, SE, P - DO NOT FILTER on MAC - USE all your data
#####
# Load BOLT and Meta data
bolt <- fread("#.bgen")
meta_wod_unfiltered <- fread("#.unfiltered")
meta_wd_unfiltered <- fread("#.unfiltered")

# Curate bolt
bolt <- bolt[, c("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF")]
setnames(bolt, c("ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF"), c("NEA", "EA", "EAF", "INFO", "BETA.BOLT", "SE.BOLT", "P.BOLT"))

# Curate wod_meta
names(meta_wod_unfiltered) <- gsub(x = names(meta_wod_unfiltered), pattern = "cohort_", replacement = "K.META_WOD")  
names(meta_wod_unfiltered) <- toupper(names(meta_wod_unfiltered))
meta_wod_unfiltered <- meta_wod_unfiltered[, c(2,1,3,5,4,7,8,9,6,10:47)]
#meta_wod_unfiltered[, !grepl("INFO", names(meta_wod_unfiltered)), with=F]
setnames(meta_wod_unfiltered, c("RSID", "POS", "ALLELE_A", "ALLELE_B", "CODED_AF", "P_VALUE", "Q", "P_HETEROGENEITY", "BETA", "SE", "I2"), 
                              c("SNP", "BP", "NEA.META_WOD", "EA.META_WOD", "EAF.META_WOD", "P.META_WOD", "Q.META_WOD", "P_Q.META_WOD", "BETA.META_WOD", "SE.META_WOD", "I2.META_WOD"))

# Curate wd_meta
names(meta_wd_unfiltered) <- gsub(x = names(meta_wd_unfiltered), pattern = "cohort_", replacement = "K.META_WD")  
names(meta_wd_unfiltered) <- toupper(names(meta_wd_unfiltered))
meta_wd_unfiltered <- meta_wd_unfiltered[, c(2,1,3,5,4,7,8,9,6,10:47)]
#meta_wd_unfiltered[, !grepl("INFO", names(meta_wd_unfiltered)), with=F]
setnames(meta_wd_unfiltered, c("RSID", "POS", "ALLELE_A", "ALLELE_B", "CODED_AF", "P_VALUE", "Q", "P_HETEROGENEITY", "BETA", "SE", "I2"), 
                              c("SNP", "BP", "NEA.META_WD", "EA.META_WD", "EAF.META_WD", "P.META_WD", "Q.META_WD", "P_Q.META_WD", "BETA.META_WD", "SE.META_WD", "I2.META_WD"))
meta_wd_unfiltered <- meta_wd_unfiltered[, -c('CHR', 'BP')]

# Generate N studies and sample-size for each SNP in the metas
meta_wod_unfiltered <- meta_wod_unfiltered %>% dplyr::mutate(n_total.META_WOD = as.integer(rowSums(dplyr::select(., ends_with('_N'))))) %>%
  dplyr::mutate(n_studies.META_WOD = as.integer(rowSums(dplyr::select(., ends_with('_N')) !=0 )))
meta_wd_unfiltered <- meta_wd_unfiltered %>% dplyr::mutate(n_total.META_WD = as.integer(rowSums(dplyr::select(., ends_with('_N'))))) %>%
  dplyr::mutate(n_studies.META_WD = as.integer(rowSums(dplyr::select(., ends_with('_N')) !=0 )))

# Flip BETAs and EAFs, as SNPTEST uses the second allele as effect allele
cols_wod_beta <- names(meta_wod_unfiltered)[grepl("BETA", names(meta_wod_unfiltered))]
meta_wod_unfiltered[, (cols_wod_beta) := lapply(.SD, function(x) x*-1), .SDcols = cols_wod_beta]
meta_wod_unfiltered[, EAF.META_WOD := 1 - EAF.META_WOD]

cols_wd_beta <- names(meta_wd_unfiltered)[grepl("BETA", names(meta_wd_unfiltered))]
meta_wd_unfiltered[, (cols_wd_beta) := lapply(.SD, function(x) x*-1), .SDcols = cols_wd_beta]
meta_wd_unfiltered[, EAF.META_WD := 1 - EAF.META_WD]

# Load Chen
chen_afr <- read.xlsx("#.xlsx", sheet=4, startRow=3) %>%
  filter(Pheno=="NEU") %>% dplyr::mutate(Author="Chen 2020") %>% 
  dplyr::select(rsID, CHR, `POS.(hg19)`, Author, EA, NEA, EAF, BETA, SE, `P-value.(GC-corrected)`) %>%
  as.data.table() %>%
  setnames(c("SNP", "CHR", "BP", "Author", "EA.Author", "NEA.Author", "EAF.Author", "BETA.Author", "SE.Author", "P.Author"))
# Load reiner
reiner <- read_excel("#.xls") %>%
  dplyr::mutate(Author="Reiner 2011") %>% 
  dplyr::select(snp, chr, position, Author, allele1, allele2, freq1, effect, stderr, pvalue) %>%
  dplyr::mutate_at(vars(allele1, allele2), toupper) %>%
  as.data.table() %>%
  setnames(c("SNP", "CHR", "BP", "Author", "EA.Author", "NEA.Author", "EAF.Author", "BETA.Author", "SE.Author", "P.Author"))
# Combine and remove duplicates
prev_gwas <- rbind(reiner, chen_afr)
prev_gwas <- prev_gwas[!duplicated(SNP),]

# Merge with my data
#prev_gwas <- merge(prev_gwas, bolt, by=c("CHR", "BP"), all.x=TRUE)
prev_gwas <- merge(prev_gwas, bolt, by=c("SNP"), all.x=TRUE)
prev_gwas <- merge(prev_gwas, meta_wod_unfiltered, by=c("SNP"), all.x=TRUE)
prev_gwas <- merge(prev_gwas, meta_wd_unfiltered, by=c("SNP"), all.x=TRUE)
prev_gwas <- prev_gwas[, -c("CHR.y", "BP.y", "CHR", "BP")]
setnames(prev_gwas, c("CHR.x", "BP.x"), c("CHR", "BP (GRCh37)"))
prev_gwas$P.Author <- as.numeric(prev_gwas$P.Author)
prev_gwas <- prev_gwas %>% dplyr::mutate(across(contains("BETA") | contains("SE") | contains("INFO"), ~round(., 3)))

fwrite(prev_gwas, "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Correlation test between EAF and BETA
cor.test(prev_gwas$EAF.Author, prev_gwas$EAF, method="spearman")
cor.test(prev_gwas$BETA.Author, prev_gwas$BETA.BOLT, method="spearman")

cor_test_eaf_afr <- ggscatter(prev_gwas, x = "EAF", y = "EAF.Author",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.y = 0.02, label.sep = "\n")
   )
cor_test_beta_afr <- ggscatter(prev_gwas, x = "BETA.BOLT", y = "BETA.Author",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.y = -1, label.sep = "\n")
   )


#####
##### 2B) make a supplementary table for ALL of the SNPs that Chen and Reiner (PMID: 21738479; 1 more here at CXCL2) 
##### identifed and report if hits are novel
#####
# Load clump table with locus leads
sumstats_table <- fread("#.txt")
sumstats_table <- sumstats_table[, c("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF", "r0.6_clump_lead", "r0.1_clump_lead", "r0.001_clump_lead", "cojo_index")]
sumstats_table <- sumstats_table[sumstats_table$SNP %in% bolt[P.BOLT < 5e-8,]$SNP ,]
setnames(sumstats_table, c("ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF"), c("NEA", "EA", "EAF", "INFO", "BETA.BOLT", "SE.BOLT", "P.BOLT"))
sumstats_table[, very_rare_variant := ifelse(EAF < 0.01, "Yes", "No")]

# Keep only leading SNPs
leads_table <- sumstats_table[r0.6_clump_lead=="Yes", ]

# Merge sumstats_table with meta_wod_unfiltered
sumstats_table <- merge(sumstats_table, meta_wod_unfiltered, by="SNP")
sumstats_table <- merge(sumstats_table, meta_wd_unfiltered, by="SNP")
sumstats_table <- sumstats_table[, -c("CHR.y", "BP.y", "CHR", "BP")]
setnames(sumstats_table, c("CHR.x", "BP.x"), c("CHR", "BP (GRCh37)"))
fwrite(sumstats_table, "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Map lead SNPs to genes
source("#.r")
snp_data <- list()
for(x in 1:length(leads_table$SNP)) {
  tryCatch({
    snp_data[[x]] <- ncbi_snp_query(leads_table$SNP[x], version="37")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
snp_data <- as.data.table(rbindlist(snp_data))
snp_data <- snp_data[, -c("maf_population")]
leads_table <- merge(leads_table, snp_data[, c("query", "gene")], by.x="SNP", by.y="query", all.x=TRUE)
leads_table[leads_table$gene=="", ]$gene <- NA

#####
##### 3) Make another supplementary table doing the same for all SNPs associated with neutrophil count identified in UKB or any European GWAS.
#####
astle_all <- fread("#.tsv")
astle <- read.xlsx("#.xlsx", sheet=1, startRow=2) %>%
  filter(!grepl("NEUT#",`Blood.Indices.Conditionally.Significantly.Associated.with.Locus`)) %>% dplyr::mutate(Author="Astle 2016") %>% 
  dplyr::select(`rsID.of.Locus.Lead.(where.available)`, `Chr.(GRCh37)`, `BP.(GRCh37)`, Author) %>%
  as.data.table() %>%
  setnames(c("SNP", "CHR", "BP", "Author"))
astle <- merge(astle, astle_all, by.x="SNP", by.y="variant_id")
astle <- astle[, c(1:4, 9, 8, 16, 12, 13, 14)]
setnames(astle, c("effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value"), 
                c("EA.Author", "NEA.Author", "EAF.Author", "BETA.Author", "SE.Author", "P.Author"))

chen_eur <- read.xlsx("#.xlsx", sheet=2, startRow=3) %>%
  filter(Pheno=="NEU") %>% dplyr::mutate(Author="Chen 2020") %>% 
  dplyr::select(rsID, CHR, `POS.(hg19)`, Author, EA, NEA, EAF, BETA, SE, `P-value.(GC-corrected)`) %>%
  as.data.table() %>%
  setnames(c("SNP", "CHR", "BP", "Author", "EA.Author", "NEA.Author", "EAF.Author", "BETA.Author", "SE.Author", "P.Author"))

# Combine and remove duplicates
prev_gwas_eur <- rbind(astle, chen_eur)
prev_gwas_eur <- prev_gwas_eur[!duplicated(SNP),]
prev_gwas_eur$P.Author <- as.numeric(prev_gwas_eur$P.Author)

# Merge with my data
#prev_gwas <- merge(prev_gwas, bolt, by=c("CHR", "BP"), all.x=TRUE)
prev_gwas_eur <- merge(prev_gwas_eur, bolt, by=c("SNP"), all.x=TRUE)
prev_gwas_eur <- merge(prev_gwas_eur, meta_wod_unfiltered, by=c("SNP"), all.x=TRUE)
prev_gwas_eur <- merge(prev_gwas_eur, meta_wd_unfiltered, by=c("SNP"), all.x=TRUE)
prev_gwas_eur <- prev_gwas_eur[, -c("CHR.y", "BP.y", "CHR", "BP")]
setnames(prev_gwas_eur, c("CHR.x", "BP.x"), c("CHR", "BP (GRCh37)"))
prev_gwas_eur <- prev_gwas_eur %>% dplyr::mutate(across(contains("BETA") | contains("SE") | contains("INFO"), ~round(., 3)))
prev_gwas_eur <- prev_gwas_eur[!is.na(SNP),]

fwrite(prev_gwas_eur, "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Create cor.plots here too, and combine with the AFR ones from above
# Save everything to file
cor_test_eaf_eur <- ggscatter(prev_gwas_eur, x = "EAF", y = "EAF.Author",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.y = 0.02, label.sep = "\n")
   )
cor_test_beta_eur <- ggscatter(prev_gwas_eur, x = "BETA.BOLT", y = "BETA.Author",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.y = -0.50, label.sep = "\n")
   )

cor_test_all <- ggarrange(cor_test_beta_afr, cor_test_eaf_afr, cor_test_beta_eur, cor_test_eaf_eur, 
  ncol = 2, nrow = 2, common.legend = FALSE, labels=c("A", "B", "C", "D"))

ggsave(plot = cor_test_all, "#.tiff", device="tiff", width = 210, height = 210, dpi=300, compression = "lzw", units="mm")

## Combine all EUR and AFR GWAS results and remove duplicates
# Pull out all SNPs from the chen_afr and astle datasets
astle_all <- fread("#.tsv")
chen_all <- fread("#.txt")
astle_all_sig <- astle_all[p_value < 5e-8,] %>% dplyr::mutate(Author="Astle 2016") %>% 
  dplyr::select(`variant`, `chromosome`, `base_pair_location`, Author, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value) %>%
  as.data.table() %>%
  setnames(c("SNP", "CHR", "BP", "Author", "EA.Author", "NEA.Author", "EAF.Author", "BETA.Author", "SE.Author", "P.Author"))

chen_all_sig <- chen_all[`p-value` < 5e-8] %>% dplyr::mutate(Author="Chen 2020") %>% 
  dplyr::select(`ID`, `CHR`, `BP`, Author, reference_allele, other_allele, eaf, beta, se, `p-value`) %>%
  as.data.table() %>%
  setnames(c("SNP", "CHR", "BP", "Author", "EA.Author", "NEA.Author", "EAF.Author", "BETA.Author", "SE.Author", "P.Author"))

novel_gwas <- rbind(reiner, astle, chen_afr, chen_eur, astle_all_sig, chen_all_sig)
novel_gwas <- novel_gwas[!duplicated(SNP),]
novel_gwas <- merge(leads_table, novel_gwas, by=c("SNP"), all.x=TRUE)
novel_gwas <- novel_gwas[, -c("CHR.y", "BP.y", "CHR", "BP")]
setnames(novel_gwas, c("CHR.x", "BP.x"), c("CHR", "BP (GRCh37)"))
novel_gwas[, Novel := ifelse(is.na(Author), "Yes", "No")]
novel_gwas <- cSplit(novel_gwas, "gene", sep = "/", direction = "long")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
IDs <- novel_gwas$gene[!is.na(novel_gwas$gene)]
genedesc <- getBM(attributes=c('external_gene_name', 'description', 'ensembl_gene_id'), filters = 'external_gene_name', values = IDs, mart =ensembl)
genedesc <- as.data.table(genedesc)
novel_gwas <- merge(novel_gwas, genedesc, by.x="gene", by.y="external_gene_name", all.x=TRUE)
novel_gwas <- novel_gwas[, c("SNP", "CHR", "BP (GRCh37)", "NEA", "EA", "EAF", "BETA.BOLT", "SE.BOLT", "P.BOLT", "r0.6_clump_lead", "r0.1_clump_lead", "r0.001_clump_lead", "cojo_index", "Novel", "gene", "description", "ensembl_gene_id")]
novel_gwas <- novel_gwas[order(P.BOLT)]

fwrite(novel_gwas, "#.txt", quote=F, row.names=F, na="", sep="\t")


#####
##### 4) Also pull my GWAS index SNPs and look for them inside the Chen SNPs
#####
chen_all <- fread("#.txt")
chen_all <- chen_all[, c("ID", "beta", "se", "p-value")]
setnames(chen_all, c("SNP", "BETA.Chen", "SE.Chen", "P.CHEN"))
bolt_clumps <- fread("#.txt")
bolt_clumps <- bolt_clumps[r0.001_clump_lead=="Yes",]
bolt_clumps[, snpid := paste0(CHR, ":", BP, "_", ALLELE1, "_", ALLELE0)]
bolt_clumps <- merge(bolt_clumps, chen_all, by="SNP", all.x=TRUE)
bolt_clumps <- bolt_clumps[, c(1:3, 5:8,11,12,14,17:24)]
setnames(bolt_clumps, c("BP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "P_BOLT_LMM_INF"), c("BP (GRCh37)", "NEA", "EA", "EAF", "BETA.BOLT", "SE.BOLT", "P.BOLT"))

fwrite(bolt_clumps, "#.txt", quote=F, row.names=F, na="", sep="\t")


#####
##### 5) Make a third supplementary table doing the same for all INDEX (COJO) SNP associated in your BOLT or a Meta Analysis - include the COJO jma estimate
#####
cojo_all <- fread("#.txt")
cojo_all <- cojo_all[, c("SNP", "shared")]
cojo_bolt <- fread("#.cojo")
cojo_bolt <- cojo_bolt[, c(2,5,6,7,8,9,11,12,13,14)]
names(cojo_bolt)[-1] <- paste0(names(cojo_bolt)[-1], ".BOLT")

cojo_meta_wod <- fread("#.cojo")
cojo_meta_wod <- cojo_meta_wod[, c(2,5,6,7,8,9,11,12,13,14)]
names(cojo_meta_wod)[-1] <- paste0(names(cojo_meta_wod)[-1], ".META_WOD")

cojo_meta_wd <- fread("#.cojo")
cojo_meta_wd <- cojo_meta_wd[, c(2,5,6,7,8,9,11,12,13,14)]
names(cojo_meta_wd)[-1] <- paste0(names(cojo_meta_wd)[-1], ".META_WD")

cojo_all <- merge(cojo_all, cojo_bolt, by="SNP", all.x=TRUE)
cojo_all <- merge(cojo_all, cojo_meta_wod, by="SNP", all.x=TRUE)
cojo_all <- merge(cojo_all, cojo_meta_wd, by="SNP", all.x=TRUE)

cojo_all <- merge(cojo_all, meta_wod_unfiltered[, c(1,4:47)], by="SNP", all.x=TRUE)
cojo_all <- merge(cojo_all, meta_wd_unfiltered[, c(1,4:47)], by="SNP", all.x=TRUE)

cojo_all <- cojo_all[shared == "bolt|NA|NA" | shared == "bolt|NA|snptest_wod_meta",][, shared := NULL]

fwrite(cojo_all, "#.txt", quote=F, row.names=F, na="", sep="\t")

#####
##### 6) Make a fourth for your SECONDARY SNPs - report the GWAS beta, SE, and P for all k-pops, meta, BOLT, AND the COJO jma estimate.
#####
cojo_secondary <- cojo_all[p.BOLT > 5e-8 | p.META_WOD > 5e-8 | p.META_WD > 5e-8 ,]
fwrite(cojo_secondary, "#.txt", quote=F, row.names=F, na="", sep="\t")


#####
##### 7) Make a LocusZoom using your BOLT data for your BOLT INDEX SNPs
#####
# plink LD generation
bolt_clumps <- fread("#.txt")
bolt_clumps <- bolt_clumps[cojo_index=="Yes",]
clump_snps <- bolt_clumps$SNP
writeLines(clump_snps, "#.txt")

bed_path <- paste0("#")
tmp <- tempfile()
system(glue::glue("plink --bfile {bed_path}",
                " --snps-only",
                " --r2",
                " --ld-window 9999999",
                " --ld-window-kb 500",
                " --ld-snp-list #.txt",
                " --ld-window-r2 0",
                " --threads 8",
                " --out {tmp}"))
ld_plink <- data.table::fread(paste0(tmp, ".ld"), data.table = TRUE, fill=TRUE)

# Load BOLT GWAS results
bolt <- fread("#.bgen")
bolt <- bolt[, c("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF")]
setnames(bolt, c("ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF"), c("EA", "NEA", "EAF", "INFO", "BETA", "SE", "P"))
bolt <- bolt[order(bolt$SNP, abs(bolt$P) ), ]
bolt <- bolt[ !duplicated(bolt$SNP), ]

snps_of_interest <- unique(c(unique(ld_plink$SNP_A), unique(ld_plink$SNP_B)))
bolt_lz <- bolt[SNP %in% snps_of_interest,]

# Load unique genes
UCSC_GRCh37_Genes_UniqueList.txt <- fread("#.txt", stringsAsFactors = FALSE, header = TRUE)

# Load the locuszoom function into R
source("#.R")

clump_snps <- clump_snps[!clump_snps %in% "8:1286073_GC_G"]
for(i in 1:length(clump_snps)) {
	locus.zoom(data = bolt_lz,
			   snp = clump_snps[i],
			   ld.file = ld_plink,
			   offset_bp = 100000,
			   genes.data = UCSC_GRCh37_Genes_UniqueList.txt,
			   plot.title = clump_snps[i],
			   file.name = paste0("#/YY_LZ_", clump_snps[i], ".tiff"),
         plot.type = "tiff",
			   secondary.label = FALSE,
			   ignore.lead = TRUE)
}

#####
##### 8) Make a publication read forrest plot
##### - 8a)  make one with meta_wd and meta_wod as a sensitivity analysis for SupMats
##### - 8b) make one excluding meta_wd for the main text.
#####
cojo_index <- fread("#.txt")
cojo_index <- data.table::melt(cojo_index, id.vars = c("SNP"), measure = list(c("b.BOLT","BETA.META_WOD","BETA.META_WD"), 
  c("se.BOLT","SE.META_WOD","SE.META_WD"), 
  c("p.BOLT","P.META_WOD","P.META_WD")), 
  value.name = c("BETA", "SE", "P-value"),
  variable.name = c("GWAS"))
cojo_index[, GWAS := fcase(
  GWAS == 1, "BOLT-LMM",
  GWAS == 2, "META-WOD",
  GWAS == 3, "META-WD"
)][, LowerLimit := BETA-1.96*SE][, UpperLimit := BETA+1.96*SE]

p = ggplot(data=cojo_index,
    aes(x = GWAS,y = BETA, ymin = LowerLimit, ymax = UpperLimit ))+
    geom_pointrange(aes(col=GWAS))+
    geom_hline(aes(fill=GWAS),yintercept =0, linetype=2)+
    xlab('Variant ID')+ ylab(expression("Variant "*beta*" coef. (95% CI)"))+
    geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=GWAS),width=0.5,cex=0.2) + expand_limits(y = c(-0.55, 1.1)) +
    facet_wrap(~SNP,strip.position="top",nrow=9,scales = "free_y") +
    theme_pubr() + scale_x_discrete(limits=rev) +
    theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    coord_flip()
ggsave(plot = p, width=180, height=210, units="mm", compression = "lzw", device = "tiff",
            dpi=300, filename = "#/YZ-forest_sup.tiff")


p = ggplot(data=cojo_index[GWAS != "META-WD",],
    aes(x = GWAS,y = BETA, ymin = LowerLimit, ymax = UpperLimit ))+
    geom_pointrange(aes(col=GWAS))+
    geom_hline(aes(fill=GWAS),yintercept =0, linetype=2)+
    xlab('Variant ID')+ ylab(expression("Variant "*beta*" coef. (95% CI)"))+
    geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=GWAS),width=0.5,cex=0.2) + expand_limits(y = c(-0.55, 1.1)) +
    facet_wrap(~SNP,strip.position="top",nrow=9,scales = "free_y") +
    theme_pubr() + scale_x_discrete(limits=rev) +
    theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    coord_flip()
ggsave(plot = p, width=180, height=210, units="mm", compression = "lzw", device = "tiff",
            dpi=300, filename = "#.tiff")


#####
##### 9) Make a Main text Table with your BOLT index SNPs and BOLT estimates, include the meta estimates where possible indicating the number of k-pops and N in the meta
#####
# Go up and load meta_wod_unfiltered and meta_wd_unfiltered
cojo_index <- fread("#.txt")
cojo_index <- cojo_index[, c(1,3,4,5,32,74,115,116)]
cojo_index <- cojo_index[order(p.BOLT)]
names(cojo_index) <- c("SNP", "BETA-BOLT", "SE-BOLT", "P-BOLT", "BETA-META-WOD", "BETA-META-WD", "META-N", "META-N-Studies")

cojo_index <- cojo_index %>% dplyr::mutate_at(vars(`BETA-BOLT`, `SE-BOLT`, `BETA-META-WOD`, `BETA-META-WD`), list(~ round(., 2)))
cojo_index[, `P-BOLT` := format.pval(`P-BOLT`, digits=3, eps=2.220446e-256)]

cojo_index_tb <- flextable(cojo_index) %>% 
    autofit() %>% 
    FitFlextableToPage()

save_as_docx(cojo_index_tb, path = "#.docx")

#####
##### 10) Make a manhattan plot with BOLT results
##### -10a) Just BOLT, to add in main text
##### -10b) BOLT vs Chen, compare AFRvAFR
##### -10c) BOLT vs Astle, compare AFRvEUR
##### Go up and load astle eur, chen afr and bolt
# Load BOLT and Meta data
bolt <- fread("#.bgen")

# Curate bolt
bolt <- bolt[, c("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF")]
setnames(bolt, c("ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF"), c("EA", "NEA", "EAF", "INFO", "BETA.BOLT", "SE.BOLT", "P.BOLT"))

astle_all <- fread("#.tsv")
chen_all <- fread("#.txt")
astle_all <- astle_all[, c("variant_id", "chromosome", "base_pair_location", "p_value")]
names(astle_all) <- c("SNP", "CHR", "POS", "pvalue")
chen_all <- chen_all[, c("ID", "CHR", "BP", "p-value")]
names(chen_all) <- c("SNP", "CHR", "POS", "pvalue")
setnames(bolt, c("BP", "P.BOLT"), c("POS", "pvalue"))
bolt <- bolt[, -c(4:9)]
plink_bim <- fread("#")
bolt <- bolt[SNP %in% plink_bim$V2, ]

bolt_clumps <- fread("#.txt")
bolt_clumps <- bolt_clumps[cojo_index=="Yes",]
clump_snps <- bolt_clumps$SNP

prep_manhattan <- function(gwas_data_load) {
  sig_data <- gwas_data_load %>% 
    subset(`pvalue` < 0.05)
  notsig_data <- gwas_data_load %>% 
    subset(`pvalue` >= 0.05) %>%
    group_by(CHR) %>% 
    sample_frac(0.1)
  gwas_data <- as.data.table(bind_rows(sig_data, notsig_data))
  return(gwas_data)
}

bolt_p <- prep_manhattan(bolt)
astle_all_p <- prep_manhattan(astle_all)
chen_all_p <- prep_manhattan(chen_all)
chen_all_p[pvalue == 0, pvalue := 1.67E-107]

## qqman
bolt_p_qqman <- bolt_p %>% setnames(c("SNP", "CHR", "BP", "P"))
astle_all_p_qqman <-astle_all_p %>% setnames(c("SNP", "CHR", "BP", "P"))
chen_all_p_qqman <- chen_all_p %>% setnames(c("SNP", "CHR", "BP", "P"))

# bolt
tiff(type="cairo", width=297, height=95, units="mm", compression = "lzw", 
            res=300, filename = "#.tiff")
manhattan(bolt_p_qqman, suggestiveline=FALSE, highlight=clump_snps, cex.axis=0.8,las=2,font=4)
dev.off()

# bolt v chen
tiff(type="cairo", width=297, height=190, units="mm", compression = "lzw", 
            res=300, filename = "#.tiff")
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(bolt_p_qqman,ylim=c(0,100),cex.axis=0.8,las=2,font=4)
par(mar=c(5,5,3,3))
manhattan(chen_all_p_qqman,ylim=c(140,0),cex.axis=0.8,las=2,font=4,xlab="",xaxt="n")
dev.off()

# bolt v chen, zoomed
tiff(type="cairo", width=297, height=210, units="mm", compression = "lzw", 
            res=300, filename = "#.tiff")
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(bolt_p_qqman,ylim=c(0,20),cex.axis=0.8,las=2,font=4)
par(mar=c(5,5,3,3))
manhattan(chen_all_p_qqman,ylim=c(20,0),cex.axis=0.8,las=2,font=4,xlab="",xaxt="n")
dev.off()

# bolt v astle
tiff(type="cairo", width=297, height=210, units="mm", compression = "lzw", 
            res=300, filename = "#.tiff")
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(bolt_p_qqman,ylim=c(0,100),cex.axis=0.8,las=2,font=4)
par(mar=c(5,5,3,3))
manhattan(astle_all_p_qqman,ylim=c(240,0),cex.axis=0.8,las=2,font=4,xlab="",xaxt="n")
dev.off()

# bolt v astle, zoomed
tiff(type="cairo", width=297, height=210, units="mm", compression = "lzw", 
            res=300, filename = "#.tiff")
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(bolt_p_qqman,ylim=c(0,20),cex.axis=0.8,las=2,font=4)
par(mar=c(5,5,3,3))
manhattan(astle_all_p_qqman,ylim=c(20,0),cex.axis=0.8,las=2,font=4,xlab="",xaxt="n")
dev.off()


#####
##### 11) Extract ALL your study index SNPs into a dosage file for possible followup sensitivity analysis in R, for example using fewer PCs, or estimating EAF by place of birth.
##### Use script from prep_guillaume
# Get dosage data using PLINK2
bolt_clumps <- fread("#.txt")
bolt_clumps <- bolt_clumps[r0.001_clump_lead=="Yes",]
clump_snps <- bolt_clumps$SNP
clump_snps <- c(clump_snps, "rs334")
writeLines(clump_snps, "#.txt")

dos_list <- list()
for (CHR in c(1,2,5,6,9,11,12,16)) {
    x <- CHR
    if(CHR < 10) CHR=paste0("0",CHR)
    bgen_file <- paste0("#.chr", CHR, ".bgen")
    sample_file <- paste0("#.chr", CHR, ".sample")
    tmp <- tempfile()
    system(glue::glue("plink2 --bgen {bgen_file} --sample {sample_file}",
                    " --extract #.txt",
                    " --export A 'include-alt'",
                    " --keep #.txt",
                    " --threads 8",
                    " --out {tmp}"))
    dos <- data.table::fread(paste0(tmp, ".raw"), data.table = TRUE)
    dos <- dos[, -c("PAT", "MAT", "SEX", "PHENOTYPE")]
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

#####
##### 12) Sensitivity BOLT GWAS
##### 
bolt <- fread("#.bgen")

bolt_sens <- fread("#.bgen.gz")
bolt_sens_sig <- bolt_sens[P_BOLT_LMM_INF < 5e-8,]

bolt_clump_cojo <- fread("#.txt")
bolt_clump_cojo <- merge(bolt_clump_cojo, bolt_sens[, c("SNP", "BETA", "SE", "P_BOLT_LMM_INF")], by="SNP", all.x=TRUE, suffixes=c("BOLT", "BOLT_sensitivity"))
setnames(bolt_clump_cojo, c("BETA", "SE", "P_BOLT_LMM_INF"), c("BETA.BOLT.sensitivity", "SE.BOLT.sensitivity", "P.BOLT.sensitivity"))

fwrite(bolt_clump_cojo, "#.txt", quote=F, row.names=F, na=NA, sep="\t")


# Write bolt sens sig
bolt_sens_sig <- bolt_sens_sig[, c("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF")]
setnames(bolt_sens_sig, c("ALLELE1", "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P_BOLT_LMM_INF"), c("NEA", "EA", "EAF", "INFO", "BETA.BOLT.sensitivity", "SE.BOLT.sensitivity", "P.BOLT.sensitivity"))

fwrite(bolt_sens_sig, "#.txt", quote=F, row.names=F, na=NA, sep="\t")





