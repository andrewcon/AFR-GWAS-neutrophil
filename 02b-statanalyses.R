# Load packages------
load_pkg <- rlang::quos(tidyverse, data.table, ggplot2, ggpubr, ggrepel, RColorBrewer, car, speedglm, foreach, doParallel, kableExtra, meta)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

# Load data------
covars <- fread("#.txt", stringsAsFactors = F)
covars <- covars[, c("genid", "neutrophil_count", "nc_log", "neutrophil_device", 
                     "sample_year", "sample_month", "sample_day", "sample_day_minutes", 
                     "sex", "age", "menopause", "cob_un_regions", "Kpop", 
                     paste0("PC", 1:100), "assessment_center", "BMI", "height", "weight", "home_north", 
                     "home_east", "townsendDI", "smoking_status", "drinker_status", 
                     "rs2814778_T")]
covars <- covars[!is.na(neutrophil_count), ]
for (i in names(covars)) {
  covars[get(i) == " " | get(i) == "", (i) := NA]
}
covars$menopause <- as.character(covars$menopause)
covars[is.na(menopause), menopause := "male"]

# Boxplots of neutrophil count by Kpop
kpop_boxplot <- ggboxplot(covars[order(Kpop)], "Kpop", "nc_log",
     fill = "Kpop", ylab = "log(Neutrophil count)")
all_histogram <- gghistogram(covars, x="nc_log", fill = "lightgray",
  xlab="Count", ylab = "log(Neutrophil count)", labels=c("A", 'B')) + 
  annotate("text", x = -2, y=rev(c(600, 700, 800, 900, 1000)), label = c(paste0("N = ",nrow(covars)), 
                                                               paste0("SD = ", signif(sd(covars$nc_log), 2) ),
                                                               paste0("Mean = ", signif(mean(covars$nc_log), 2) ),
                                                               paste0("Median = ", signif(median(covars$nc_log), 2) ),
                                                               paste0("Units = log(", expression(10^9), " cells/Litre)") ),
                                                               hjust = 0)
nc_plot <- ggarrange(all_histogram, kpop_boxplot, ncol=1, nrow=2, labels=c("A", "B"))
ggsave(plot = nc_plot, 
       filename = "#.tiff", 
        device="tiff", width = 210, height = 210, dpi=300, compression = "lzw", units="mm")

# Anova type I, II and III------
# Remove NA from covars
covars <- na.omit(covars)
nr_na <- 5976 - 5706
#covars <- covars[!rowSums(covars == "prefer_not_answer") > 0, ]
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
#covars$rm_astle <- factor(covars$rm_astle, 
#                          levels = c(0, 1), 
#                          labels = c("No", "Yes")
#)
covars <- covars[, -c("rm_astle")]
covars$smoking_status <- factor(covars$smoking_status, 
                                levels = c("prefer_not_answer", "never", "previous", "current"), 
                                labels = c("Prefer not answer", "Never", "Previous", "Current")
)
covars$drinker_status <- factor(covars$drinker_status, 
                                levels = c("prefer_not_answer", "never", "previous", "current"), 
                                labels = c("Prefer not answer", "Never", "Previous", "Current")
)
covars$rs2814778_T <- factor(covars$rs2814778_T, 
                             levels = c(0, 1, 2), 
                             labels = c("TT", "TC", "CC")
)
covars$assessment_center <- as.factor(covars$assessment_center)
lm_formula <- colnames(covars)[-c(1, 2, 3, 24:113)]
lm_formula <- paste0(lm_formula, collapse = " + ")
lm_formula <- as.formula(paste0("nc_log ~ ", lm_formula))
covars_lm <- lm(lm_formula, data = covars)
# Univariate analysis
covars.univariate <- data.table()
x <- colnames(covars)[-c(1, 2, 3, 24:113)]
for(v in x) {
  temp_formula <- as.formula(paste0("nc_log ~ ", v))
  temp_lm <- lm(temp_formula, data = covars)
  temp_aov <- aov(temp_lm)
  temp_summary <- summary(temp_aov)
  afss <- sum(temp_summary[[1]][[2]])
  temp.univariate <- data.table(
    "rn" = v, 
    "DF" = temp_summary[[1]][1,1],
    "Sum Sq" = temp_summary[[1]][1,2],
    "Mean Sq" = temp_summary[[1]][1,3],
    "F value" = temp_summary[[1]][1,4],
    "PctExp" = temp_summary[[1]][1,2]/afss*100,
    "type" = "univariate"
  )
  covars.univariate <- rbind(covars.univariate, temp.univariate)
}
# Type I ANOVA
covars.I.aov <- aov(covars_lm)
covars.I.aov <- summary(covars.I.aov)[[1]]
covars.I.aov <- setDT(covars.I.aov, keep.rownames = TRUE)[]
afss <- covars.I.aov$"Sum Sq"
covars.I.aov[, PctExp := afss/sum(afss)*100]
covars.I.aov$type <- "anova_one"
covars.I.aov[, "rn" := lapply(.SD, trimws), .SDcols = "rn"]
# Type II ANOVA
covars.II.aov <- car::Anova(covars_lm, type = 2)
covars.II.aov <- setDT(covars.II.aov, keep.rownames = TRUE)[]
afss <- covars.II.aov$"Sum Sq"
covars.II.aov[, PctExp := afss/sum(afss)*100]
covars.II.aov$type <- "anova_two"
# Merge results
covars.all.aov <- rbind(covars.univariate, covars.I.aov, covars.II.aov, fill = TRUE)
covars.all.aov$PctExp <- signif(covars.all.aov$PctExp, digits = 3)
covars.all.aov <- covars.all.aov[!startsWith(rn, "Residuals")]
covars.all.aov <- covars.all.aov[dim(covars.all.aov)[1]:1,]
x <- colnames(covars)[-c(1, 2, 3, 24:113)]
covars.all.aov$type <- factor(covars.all.aov$type, 
                              levels = c("anova_two", "anova_one", "univariate"),
                              labels = c("ANOVA II", "ANOVA I (top-bottom order)", "Univariate"))
# Plot variance explained
# Define scale
aov_scale_manual <- scale_fill_manual(
  "Analysis type",
  values = c("#4DAF4A", "#377EB8", "#E41A1C"),
  guide = guide_legend(reverse = TRUE)
)
aov_scale_x <- scale_x_discrete(
  breaks = x,
  labels = c("Neutrophil Device", "Sample year", "Sample month", "Sample day", "Sample day\nminutes",
             "Sex", "Age", "Menopause", "Continents &\nWorld UN regions", "K-mean cluster",
             "PC1", "PC2", "PC3", "PC4", "PC5",
             "PC6", "PC7", "PC8", "PC9", "PC10", "Assessment Centre",
             "Body mass index", "Height", "Weight", "Home north coord.",
             "Home east coord.", "Depravation index", "Smoking status", "Drinker status",
             "rs2814778 C allele")
)
v_barplot <- ggbarplot(
  data = covars.all.aov,
  x = "rn",
  y = "PctExp",
  combine = TRUE,
  color = "black",
  fill = "type",
  palette = NULL,
  orientation = "horiz",
  size = 0.25,
  width = 0.8,
  title = NULL,
  xlab = "",
  ylab = "Variance explained %",
  label = TRUE,
  lab.col = "black",
  lab.size = 3,
  lab.pos = c("out"),
  lab.vjust = 0.45,
  lab.hjust = -0.2,
  lab.nb.digits = 2,
  sort.by.groups = TRUE,
  position = position_dodge(0.9)
) + aov_scale_manual + aov_scale_x +
  theme(axis.text.y = element_text(size=10),
        legend.position = c(0.80, 0.925),
        legend.direction = "vertical",
        legend.background = element_rect(colour = "black", size = 0.5)) + 
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) + 
  guides(colour = guide_legend(reverse=TRUE))

ggsave(plot = v_barplot, 
       filename = "#.tiff", 
       device="tiff", width = 210, height = 290, dpi=300, compression = "lzw", units="mm")

##########################################################
# View how Duffy allele affects neutrophil count by Kpop #
##########################################################
###
### 1) Without adjusting for PCs
###
lm_formula_noPC <- c("sex", "neutrophil_device", "sample_year", "sample_month",
  "sample_day", "sample_day_minutes", "age", "assessment_center", "rs2814778_T")
lm_formula_noPC <- paste0(lm_formula_noPC, collapse = " + ")
lm_formula_noPC <- formula(paste0("nc_log ~ ", lm_formula_noPC))

kpop_info <- data.table(Kpop = character(), sample_size = numeric(), analysis_type = character(), nc_log_mean = numeric(), nc_log_median = numeric(), nc_log_sd = numeric(), 
  beta_rs2814778_TC = numeric(), se_rs2814778_TC = numeric(), pval_rs2814778_TC = numeric(), count_rs2814778_TC = numeric(),
  beta_rs2814778_CC = numeric(), se_rs2814778_CC = numeric(), pval_rs2814778_CC = numeric(), count_rs2814778_CC = numeric(),
  rs2814778_varExp = numeric())
kpops <- c("K1", "K2", "K3", "K4", "K5", "K6", "K7")

for(x in 1:7) {
  kpop_lm <- lm(lm_formula_noPC, data = covars[Kpop == kpops[x]])
  kpop_aov <- aov(kpop_lm)
  kpop_summary <- summary(kpop_aov)
  afss <- sum(kpop_summary[[1]][[2]][c(9,10)])
  kpop_dt <- data.table(summary(kpop_lm)[[4]], keep.rownames=T)
  if (length(table(covars[Kpop == kpops[x], 'rs2814778_T'])) == 3) {
    count_rs2814778_TC <- 2*table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TT"]] + table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TC"]]
    count_rs2814778_CC <- 2*table(covars[Kpop == kpops[x], 'rs2814778_T'])[["CC"]] + table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TC"]]
  } else if (length(table(covars[Kpop == kpops[x], 'rs2814778_T'])) == 2) {
    count_rs2814778_TC <- table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TC"]]
    count_rs2814778_CC <- 2*table(covars[Kpop == kpops[x], 'rs2814778_T'])[["CC"]] + table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TC"]]
  } else {
    count_rs2814778_TC <- 0
    count_rs2814778_CC <- 2*table(covars[Kpop == kpops[x], 'rs2814778_T'])[["CC"]]
  }
  kpop_row <- data.table(Kpop = kpops[x], sample_size = table(covars$Kpop)[[x]], analysis_type = "No PC adjustment", nc_log_mean = mean(covars[Kpop == kpops[x],nc_log]), nc_log_median = median(covars[Kpop == kpops[x],nc_log]), 
  nc_log_sd = sd(covars[Kpop == kpops[x], nc_log]), 
  beta_rs2814778_TC = kpop_dt[rn == "rs2814778_TTC"]$Estimate, 
  se_rs2814778_TC = kpop_dt[rn == "rs2814778_TTC"]$'Std. Error', 
  pval_rs2814778_TC = kpop_dt[rn == "rs2814778_TTC"]$'Pr(>|t|)',
  count_rs2814778_TC = count_rs2814778_TC,
  beta_rs2814778_CC = kpop_dt[rn == "rs2814778_TCC"]$Estimate, 
  se_rs2814778_CC = kpop_dt[rn == "rs2814778_TCC"]$'Std. Error', 
  pval_rs2814778_CC = kpop_dt[rn == "rs2814778_TCC"]$'Pr(>|t|)', 
  count_rs2814778_CC = count_rs2814778_CC,
  rs2814778_varExp = kpop_summary[[1]][9,2]/afss*100)
  kpop_info <- rbind(kpop_info, kpop_row)
}

###
### 2) Adjusting for PCs
###
lm_formula_noPC <- c("sex", "neutrophil_device", "sample_year", "sample_month",
  "sample_day", "sample_day_minutes", "age", "assessment_center", paste0("PC", 1:100), "rs2814778_T")
lm_formula_noPC <- paste0(lm_formula_noPC, collapse = " + ")
lm_formula_noPC <- formula(paste0("nc_log ~ ", lm_formula_noPC))

for(x in 1:7) {
  kpop_lm <- lm(lm_formula_noPC, data = covars[Kpop == kpops[x]])
  kpop_aov <- aov(kpop_lm)
  kpop_summary <- summary(kpop_aov)
  afss <- sum(kpop_summary[[1]][[2]][c(109,110)])
  kpop_dt <- data.table(summary(kpop_lm)[[4]], keep.rownames=T)
  if (length(table(covars[Kpop == kpops[x], 'rs2814778_T'])) == 3) {
    count_rs2814778_TC <- 2*table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TT"]] + table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TC"]]
    count_rs2814778_CC <- 2*table(covars[Kpop == kpops[x], 'rs2814778_T'])[["CC"]] + table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TC"]]
  } else if (length(table(covars[Kpop == kpops[x], 'rs2814778_T'])) == 2) {
    count_rs2814778_TC <- table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TC"]]
    count_rs2814778_CC <- 2*table(covars[Kpop == kpops[x], 'rs2814778_T'])[["CC"]] + table(covars[Kpop == kpops[x], 'rs2814778_T'])[["TC"]]
  } else {
    count_rs2814778_TC <- 0
    count_rs2814778_CC <- 2*table(covars[Kpop == kpops[x], 'rs2814778_T'])[["CC"]]
  }
  kpop_row <- data.table(Kpop = kpops[x], sample_size = table(covars$Kpop)[[x]], analysis_type = "With PC1:100 adjustment", nc_log_mean = mean(covars[Kpop == kpops[x],nc_log]), nc_log_median = median(covars[Kpop == kpops[x],nc_log]), 
  nc_log_sd = sd(covars[Kpop == kpops[x], nc_log]), 
  beta_rs2814778_TC = kpop_dt[rn == "rs2814778_TTC"]$Estimate, 
  se_rs2814778_TC = kpop_dt[rn == "rs2814778_TTC"]$'Std. Error', 
  pval_rs2814778_TC = kpop_dt[rn == "rs2814778_TTC"]$'Pr(>|t|)',
  count_rs2814778_TC = count_rs2814778_TC,
  beta_rs2814778_CC = kpop_dt[rn == "rs2814778_TCC"]$Estimate, 
  se_rs2814778_CC = kpop_dt[rn == "rs2814778_TCC"]$'Std. Error', 
  pval_rs2814778_CC = kpop_dt[rn == "rs2814778_TCC"]$'Pr(>|t|)', 
  count_rs2814778_CC = count_rs2814778_CC,
  rs2814778_varExp = kpop_summary[[1]][109,2]/afss*100)
  kpop_info <- rbind(kpop_info, kpop_row)
}

# Write this to observational folder in results directory
fwrite(kpop_info, "#.txt", quote=F, row.names=F, na=NA, sep="\t")


###############################################################
# Analyse Chen/Guillaume AFR SNPs associated with Neutrophil count #
###############################################################
covars_dosages <- fread("#.txt")
covars_dosages <- covars_dosages[!is.na(neutrophil_count), ]
names(covars_dosages) <- gsub(x = names(covars_dosages), pattern = "/", replacement = "_")  
names(covars_dosages) <- gsub(x = names(covars_dosages), pattern = "[()]", replacement = "")  

names(covars_dosages)[names(covars_dosages) == "rs2814778_C_T"]
covars_dosages[,"rs2814778_C_T"]

###
### 1) Whole sample, no adjustment for PCs 
###
lm_formula <- c()
cag_lm_fun <- function(x, my_lm_vars) {
  lm_formula <- my_lm_vars
  lm_formula <- paste0(lm_formula, collapse = " + ")
  lm_formula <- formula(paste0("nc_log ~ ", lm_formula))
  cag_lm <- speedlm(lm_formula, data = covars_dosages)
  cag_lm_summary <- as.data.table(summary(cag_lm)[[6]], keep.rownames=T)
  snp_rsid <- sub("_.*", "", colnames(covars_dosages)[x])
  snp_beta <- cag_lm_summary[rn == colnames(covars_dosages)[x], 2][[1]]
  snp_se <- cag_lm_summary[rn == colnames(covars_dosages)[x], 3][[1]]
  snp_pval <- cag_lm_summary[rn == colnames(covars_dosages)[x], 5][[1]]
  sumstats_snp <- data.table(pop = "CAG", rsid = snp_rsid, beta = snp_beta, se = snp_se, pval = snp_pval)
  return(sumstats_snp)
}

myCluster <- makeCluster(6, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

cag_lm_list <- list()
cag_lm_list <- foreach(snp_nr = 138:ncol(covars_dosages)) %dopar% {
  cag_lm_fun(
    x = snp_nr,
    my_lm_vars = c("sex", "neutrophil_device", "sample_year", "sample_month",
                   "sample_day", "sample_day_minutes", "age", "assessment_center", colnames(covars_dosages)[snp_nr])
  )
}
cag_lm_list <- rbindlist(cag_lm_list)

stopCluster(myCluster)

###
### 2) Whole sample, with adjustment for PCs 
###
myCluster <- makeCluster(6, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

cag_lm_list_pc <- foreach(snp_nr = 138:ncol(covars_dosages)) %dopar% {
  cag_lm_fun(
    x = snp_nr,
    my_lm_vars = c("sex", "neutrophil_device", "sample_year", "sample_month",
                   "sample_day", "sample_day_minutes", "age", "assessment_center", paste0("PC", 1:100), colnames(covars_dosages)[snp_nr])
  )
}
cag_lm_list_pc <- rbindlist(cag_lm_list_pc)

stopCluster(myCluster)

# Save summary data to file as replication
summary_chen_noPC <- summary(cag_lm_list)
summary_chen_withPC <- summary(cag_lm_list_pc)

###
### Get top 10 SNPs with lowest P-val from cag_lm_list_pc
###
top_ten_snps <- cag_lm_list_pc[order(pval),][1:10,]
top_ten_snps <- paste(top_ten_snps$rsid, collapse = "|")
covars_dosages_tts <- cbind(covars_dosages[, 1:137], covars_dosages[, grepl(top_ten_snps, names(covars_dosages)), with = F])

###
### 3) Do same analysis as in 1), but by Kpop.
### Perform meta-analysis of Kpop results
###
lm_formula <- c()
cag_lm_fun_kpop <- function(x, my_lm_vars) {
  print(paste0("x = ", x))
  print(paste0(my_lm_vars, collapse=" "))
  kpops <- c("K1", "K2", "K3", "K4", "K5", "K6", "K7")
  lm_formula <- copy(my_lm_vars)
  lm_formula <- paste0(lm_formula, collapse = " + ")
  lm_formula <- formula(paste0("nc_log ~ ", lm_formula))
  sumstats_snp_list <- list()
  for(i in 1:7){
    cag_lm <- speedlm(lm_formula, data = covars_dosages_tts[Kpop == kpops[i]])
    cag_lm_summary <- as.data.table(summary(cag_lm)[[6]], keep.rownames=T)
    snp_rsid <- sub("_.*", "", colnames(covars_dosages_tts)[x])
    snp_beta <- cag_lm_summary[rn == colnames(covars_dosages_tts)[x], 2][[1]]
    snp_se <- cag_lm_summary[rn == colnames(covars_dosages_tts)[x], 3][[1]]
    snp_pval <- cag_lm_summary[rn == colnames(covars_dosages_tts)[x], 5][[1]]
    sumstats_snp_list[[i]] <- sumstats_snp <- data.table(pop = kpops[i], rsid = snp_rsid, beta = snp_beta, se = snp_se, pval = snp_pval)
  }
  for(i in 1:7){
    covars_dosages_tts_kpop <- covars_dosages_tts[Kpop == kpops[i]]
    res_lm_formula <- copy(my_lm_vars)
    res_lm_formula <- paste0(res_lm_formula[-length(my_lm_vars)], collapse = " + ")
    res_lm_formula <- formula(paste0("nc_log ~ ", res_lm_formula))
    residuals_cag_lm <- residuals(lm(res_lm_formula, data = covars_dosages_tts_kpop))
    covars_dosages_tts_kpop$residuals <- residuals_cag_lm
    res_lm_formula <- formula(paste0("residuals ~ ", colnames(covars_dosages_tts_kpop)[x]))
    cag_lm <- speedlm(res_lm_formula, data = covars_dosages_tts_kpop)
    cag_lm_summary <- as.data.table(summary(cag_lm)[[6]], keep.rownames=T)
    snp_rsid <- sub("_.*", "", colnames(covars_dosages_tts_kpop)[x])
    snp_beta <- cag_lm_summary[rn == colnames(covars_dosages_tts_kpop)[x], 2][[1]]
    snp_se <- cag_lm_summary[rn == colnames(covars_dosages_tts_kpop)[x], 3][[1]]
    snp_pval <- cag_lm_summary[rn == colnames(covars_dosages_tts_kpop)[x], 5][[1]]
    sumstats_snp_list[[i+7]] <- sumstats_snp <- data.table(pop = paste0("residuals_", kpops[i]), rsid = snp_rsid, beta = snp_beta, se = snp_se, pval = snp_pval)
  }
  sumstats_snp <- rbindlist(sumstats_snp_list)
  return(sumstats_snp)
}

myCluster <- makeCluster(6, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

kpop_lm_list_pc <- foreach(snp_nr = 138:ncol(covars_dosages_tts)) %dopar% {
  cag_lm_fun_kpop(
    x = snp_nr,
    my_lm_vars = c("sex", "neutrophil_device", "sample_year", "sample_month",
                   "sample_day", "sample_day_minutes", "age", "assessment_center", paste0("PC", 1:100), colnames(covars_dosages_tts)[snp_nr])
  )
}
kpop_lm_list_pc <- rbindlist(kpop_lm_list_pc)

stopCluster(myCluster)

# Meta analyse by Kpop
kpop_lm_list_pc_res <- kpop_lm_list_pc[startsWith(pop, "residuals"),]
kpop_lm_list_pc <- kpop_lm_list_pc[!startsWith(pop, "residuals"),]
kpop_meta <- metagen(data = kpop_lm_list_pc, TE = kpop_lm_list_pc$beta, seTE = kpop_lm_list_pc$se, pval = kpop_lm_list_pc$pval,
                     method.tau = "REML", studlab = kpop_lm_list_pc$pop, fixed = TRUE, subgroup = kpop_lm_list_pc$rsid,
                     test.subgroup = TRUE, prediction.subgroup = TRUE, print.subgroup.name = FALSE)
kpop_meta_residuals <- metagen(data = kpop_lm_list_pc_res, TE = kpop_lm_list_pc$beta, seTE = kpop_lm_list_pc$se, pval = kpop_lm_list_pc$pval,
                     method.tau = "REML", studlab = kpop_lm_list_pc$pop, fixed = TRUE, subgroup = kpop_lm_list_pc$rsid,
                     test.subgroup = TRUE, prediction.subgroup = TRUE, print.subgroup.name = FALSE)

# Forest plot of Kpop estimates
palette_kpop <- get_palette(palette = "default", 7)
kpop_lm_duffy <- kpop_lm_list_pc[rsid == "rs2814778",]
kpop_lm_duffy[, LowerLimit := beta - se*1.96]
kpop_lm_duffy[, UpperLimit := beta + se*1.96]
p = ggplot(data=kpop_lm_duffy,
    aes(x = pop,y = beta, ymin = LowerLimit, ymax = UpperLimit ))+
    geom_pointrange(aes(col=pop))+
    geom_hline(aes(fill=pop),yintercept =0, linetype=2)+
    xlab('rs2814778')+ ylab("Beta (95% Confidence Interval)")+
    geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=pop),width=0.3,cex=1)+ 
    theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    labs(colour="K-pop") + coord_flip() + scale_x_discrete(limits=rev) +  
    scale_color_manual(values=palette_kpop) +
    scale_fill_manual(values=palette_kpop) +
    theme_pubr()
ggsave(plot = p, 
       filename = "#.tiff", 
       device="tiff", width = 210, height = 210, dpi=300, compression = "lzw", units="mm")

###
### 5) Forest plot for top 10 SNPs from dosage analysis
### Pull data from BOLT, Lettre 
###
cag_tts <- cag_lm_list_pc[order(pval),][1:10,]
kpop_meta_tts <- data.table(pop = "kpop_meta", rsid = kpop_meta[[115]], beta = kpop_meta[[116]], se = kpop_meta[[117]], pval = kpop_meta[[122]])
kpop_meta_tts_residuals <- data.table(pop = "kpop_meta_residuals", rsid = kpop_meta_residuals[[115]], beta = kpop_meta_residuals[[116]], se = kpop_meta_residuals[[117]], pval = kpop_meta_residuals[[122]])
bolt_tts <- fread(cmd=paste0("grep -wE ", "'",top_ten_snps, "'", " #.bgen"), select = c(1, 11:12, 14))
setnames(bolt_tts, c("V1", "V11", "V12", "V14"), c("rsid", "beta", "se", "pval"))
bolt_tts$pop <- "CAG_BOLT"
bolt_tts <- bolt_tts[, c(5, 1:4)]
bolt_tts$beta <- -1 * bolt_tts$beta
guillaume_tts <- fread(cmd=paste0("grep -wE ", "'",top_ten_snps, "'", " #.txt"), select = c(22, 7:9))
setnames(guillaume_tts, c("V22", "V7", "V8", "V9"), c("rsid", "beta", "se", "pval"))
guillaume_tts$pop <- "Guillaume"
guillaume_tts <- guillaume_tts[, c(5, 1:4)]
guillaume_tts$beta <- -1 * guillaume_tts$beta

tts_combined <- rbind(cag_tts, kpop_meta_tts, kpop_meta_tts_residuals, bolt_tts, guillaume_tts)
tts_combined[, LowerLimit := beta - se*1.96]
tts_combined[, UpperLimit := beta + se*1.96]
tts_combined[, pop := fcase(
  pop == "CAG", "R - Whole sample",
  pop == "kpop_meta", "R - Meta-analysis of K-pops",
  pop == "kpop_meta_residuals", "R - Meta-analysis of K-pops (residuals)",
  pop == "CAG_BOLT", "BOLT-LMM - Whole sample",
  pop == "Guillaume", "Chen et al."
)]

p <- ggplot(data=tts_combined,
    aes(x = pop,y = beta, ymin = LowerLimit, ymax = UpperLimit ))+
    geom_pointrange(aes(col=pop))+
    geom_hline(aes(fill=pop),yintercept =0, linetype=2)+
    xlab('SNP rsID')+ ylab("Beta (95% Confidence Interval)")+
    geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=pop),width=0.3,cex=0.8)+ 
    facet_wrap(~rsid,strip.position="left",nrow=9,scales = "free_y") + theme_pubr() +
    theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    labs(colour="Source") + coord_flip() + scale_x_discrete(limits=rev) + guides(colour=guide_legend(nrow=2,byrow=TRUE))
ggsave(plot = p, 
       filename = "#.tiff", 
       device="tiff", width = 210, height = 210, dpi=300, compression = "lzw", units="mm")



