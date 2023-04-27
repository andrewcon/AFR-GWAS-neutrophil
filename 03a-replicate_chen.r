# Load packages------
load_pkg <- rlang::quos(tidyverse, data.table, ggplot2, ggpubr, ggrepel, RColorBrewer, car, speedglm, foreach, doParallel, kableExtra, meta)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

###############################################################
# Analyse Chen/Guillaume AFR SNPs associated with Neutrophil count #
# So unlike in step06-d_statanalyses.R, we are using
# the same covariates that Chen et al. used
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
    my_lm_vars = c("sex", "age", "age^2", "assessment_center", colnames(covars_dosages)[snp_nr])
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
    my_lm_vars = c("sex", "age", "age^2", "assessment_center", paste0("PC", 1:10), colnames(covars_dosages)[snp_nr])
  )
}
cag_lm_list_pc <- rbindlist(cag_lm_list_pc)

stopCluster(myCluster)

# Save summary data to file as replication
summary_chen_noPC <- c(as.vector(summary(cag_lm_list$pval)))
summary_chen_noPC <- c("13139", summary_chen_noPC, cag_lm_list[pval < 0.05, .N]/13139)
summary_chen_withPC <- c(as.vector(summary(cag_lm_list_pc$pval)))
summary_chen_withPC <- c("13139", summary_chen_withPC, cag_lm_list_pc[pval < 0.05, .N]/13139)

rep_df <- data.table(
    Summary = c("N SNPs", "Min", "1st Quartile", "Median", "Mean", "3rd Quartile", "Max", "Percentage P < 0.05"),
    `P-value no PC covariates` = summary_chen_noPC,
    `P-value with PC1:10 as covariates` = summary_chen_withPC
)

fwrite(rep_df, "#.txt", na=NA, row.names=FALSE, quote=FALSE, sep="\t")

