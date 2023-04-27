#!/usr/bin/env Rscript

# Load relevant libraries------
list_of_packages <- c(
  "foreach",
  "doParallel",
  "kableExtra",
  "data.table",
  "TwoSampleMR",
  "ieugwasr",
  "tidyr",
  "rlang",
  "mr.raps"
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


# Functions------
malaria_mr <- function(exp_path, exp_name, out_path, out_name, reverse = FALSE) {
  # Load data and format it
  tryCatch({
  if(reverse != TRUE) {
    exp_dat <- fread(exp_path, stringsAsFactors = F)
    exp_dat[, BETA := exp(BETA)][, SE := exp(SE)]
    out_dat <- fread(out_path, stringsAsFactors = F)
    #Format exposure data
    exp_dat <- format_data(
      exp_dat, type="exposure", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
      eaf_col = "A1FREQ", effect_allele_col = "ALLELE0", other_allele_col = "ALLELE1", 
      pval_col = "P_BOLT_LMM_INF", chr_col = "CHR", pos_col = "BP"
    )
    # Format outcome data
    out_dat <- format_data(
      out_dat, type="outcome", snp_col = "MarkerName", beta_col = "Effect", se_col = "StdErr",
      eaf_col = "Freq1", effect_allele_col = "Allele2", other_allele_col = "Allele1", 
      pval_col = "P-value", chr_col = "Chromosome", pos_col = "Position", 
      snps = exp_dat$SNP
    )

  # Harmonize data
  dat <- harmonise_data(
    exposure_dat = exp_dat, 
    outcome_dat = out_dat
  )
  dat_lowT <- copy(dat)

  # Perform clumping for 5e-8
  setnames(dat, c("SNP", "pval.exposure"), c("rsid", "pval"))
  dat <- ld_clump(dat=dat, clump_kb=10000, clump_r2=0.001, clump_p=5e-8, 
                    pop = "AFR", access_token=NULL, bfile="#", 
                    plink_bin="#")
  setnames(dat, c("rsid", "pval"), c("SNP", "pval.exposure"))

  dat$id.exposure <- exp_name
  dat$id.outcome <- out_name
  dat <- dat[dat$SNP != "rs72713693", ]

  # Perform clumping for 5e-5
  setnames(dat_lowT, c("SNP", "pval.exposure"), c("rsid", "pval"))
  dat_lowT <- ld_clump(dat=dat_lowT, clump_kb=10000, clump_r2=0.001, clump_p=5e-5, 
                    pop = "AFR", access_token=NULL, bfile="#", 
                    plink_bin="#")
  setnames(dat_lowT, c("rsid", "pval"), c("SNP", "pval.exposure"))

  dat_lowT$id.exposure <- exp_name
  dat_lowT$id.outcome <- out_name

  # Perform MR
  res <- mr(dat)
  res_raps <- mr_raps_andrei(dat_lowT[dat_lowT$mr_keep == TRUE,])
  res <- rbind(res, res_raps)
  res <- generate_odds_ratios(res)
  res$exposure <- exp_name
  res$outcome <- out_name
  if(res$nsnp != 1) {
      het <- mr_heterogeneity(dat)
      het$exposure <- exp_name
      het$outcome <- out_name
      ple <- mr_pleiotropy_test(dat)
      ple$exposure <- exp_name
      ple$outcome <- out_name
      res_single <- mr_singlesnp(dat)
      res_single$exposure <- exp_name
      res_single$outcome <- out_name
      res_loo <- mr_leaveoneout(dat)
      res_loo$exposure <- exp_name
      res_loo$outcome <- out_name
      #dir_test <- directionality_test(dat)
      #dir_test$exposure <- exp_name
      #dir_test$outcome <- out_name
    } else {
      het <- NA
      ple <- NA
      res_single <- NA
      res_loo <- NA
      dir_test <- NA
    }
  res_all <- list(res, het, ple, res_single, res_loo, dat_lowT)
  print(paste0("Finished ", exp_name, " on ", out_name, "."))
  return(res_all)
  
  } else {
    exp_dat <- fread(exp_path, stringsAsFactors = F)
    exp_dat <- exp_dat[`P-value` < 5e-5,]
    out_dat <- fread(out_path, stringsAsFactors = F)
    out_dat[, BETA := exp(BETA)][, SE := exp(SE)]

    # Perform clumping for 5e-8
    setnames(exp_dat, c("MarkerName", "P-value"), c("rsid", "pval"))
    exp_dat <- ld_clump(dat=exp_dat, clump_kb=10000, clump_r2=0.001, clump_p=5e-8, 
                      pop = "AFR")
    setnames(exp_dat, c("rsid", "pval"), c("MarkerName", "P-value"))

    exp_dat <- format_data(
      exp_dat, type="exposure", snp_col = "MarkerName", beta_col = "Effect", se_col = "StdErr",
      eaf_col = "Freq1", effect_allele_col = "Allele2", other_allele_col = "Allele1", 
      pval_col = "P-value", chr_col = "Chromosome", pos_col = "Position", 
      snps = exp_dat$SNP
    )
    # Format outcome data
    out_dat <- format_data(
      out_dat, type="outcome", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
      eaf_col = "A1FREQ", effect_allele_col = "ALLELE0", other_allele_col = "ALLELE1", 
      pval_col = "P_BOLT_LMM_INF", chr_col = "CHR", pos_col = "BP", snps = exp_dat$SNP)
    
    
    # Harmonize data
    dat <- harmonise_data(
      exposure_dat = exp_dat, 
      outcome_dat = out_dat
    )
    dat_lowT <- copy(dat)

    dat$id.exposure <- exp_name
    dat$id.outcome <- out_name

    # Perform clumping for 5e-5
    setnames(dat_lowT, c("SNP", "pval.exposure"), c("rsid", "pval"))
    dat_lowT <- ld_clump(dat=dat_lowT, clump_kb=10000, clump_r2=0.001, clump_p=5e-5, 
                      pop = "AFR", access_token=NULL, bfile="#", 
                      plink_bin="#")
    setnames(dat_lowT, c("rsid", "pval"), c("SNP", "pval.exposure"))

    dat_lowT$id.exposure <- exp_name
    dat_lowT$id.outcome <- out_name

    # Perform MR
    res <- mr(dat)
    res_raps <- mr_raps_andrei(dat_lowT[dat_lowT$mr_keep == TRUE,])
    res <- rbind(res, res_raps)
    res <- generate_odds_ratios(res)
    res$exposure <- exp_name
    res$outcome <- out_name
    if(res$nsnp != 1) {
        het <- mr_heterogeneity(dat)
        het$exposure <- exp_name
        het$outcome <- out_name
        ple <- mr_pleiotropy_test(dat)
        ple$exposure <- exp_name
        ple$outcome <- out_name
        res_single <- mr_singlesnp(dat)
        res_single$exposure <- exp_name
        res_single$outcome <- out_name
        res_loo <- mr_leaveoneout(dat)
        res_loo$exposure <- exp_name
        res_loo$outcome <- out_name
        #dir_test <- directionality_test(dat)
        #dir_test$exposure <- exp_name
        #dir_test$outcome <- out_name
      } else {
        het <- NA
        ple <- NA
        res_single <- NA
        res_loo <- NA
        dir_test <- NA
      }
    res_all <- list(res, het, ple, res_single, res_loo, dat)
    print(paste0("Finished ", exp_name, " on ", out_name, "."))
    return(res_all)
  }
  }, error=function(e){})
}

mr_raps_andrei <- function(mydat, over.dispersion=FALSE, loss.function="huber") {
    data <- data.frame(beta.exposure = mydat$beta.exposure,
                       beta.outcome = mydat$beta.outcome,
                       se.exposure = mydat$se.exposure,
                       se.outcome = mydat$se.outcome)
    out <- suppressMessages(
        mr.raps(b_exp = mydat$beta.exposure, 
                         b_out = mydat$beta.outcome, 
                         se_exp = mydat$se.exposure, 
                         se_out = mydat$se.outcome,
                         over.dispersion = over.dispersion,
                         loss.function = loss.function))
    res_list <- list(b = out$beta.hat,
         se = out$beta.se,
         pval = pnorm(- abs(out$beta.hat / out$beta.se)) * 2,
         nsnp = length(mydat$beta.exposure))
    mr_tab <- data.frame(
    	id.exposure = mydat$id.exposure[1],
			id.outcome = mydat$id.outcome[1],
      outcome = "outcome",
      exposure = "exposure",
			method = "MR RAPS",
			nsnp = res_list[[4]],
			b = res_list[[1]],
			se = res_list[[2]],
			pval = res_list[[3]]
		)
		return(mr_tab)
}

# Load the neutrophil and malaria datasets------
# Create table with MR information
mr_info <- data.table(
    exposure_name = c("Neutrophil", "Neutrophil", "Neutrophil", "Neutrophil", "Malaria_Overall" , "Malaria_CM", "Malaria_SMA", "Malaria_Other"),
    exposure_path = c(
        "#.txt", 
        "#.txt", 
        "#.txt", 
        "#.txt", 
        "#.txt", 
        "#.txt", 
        "#.txt", 
        "#.txt"),
    outcome_name = c("Malaria_Overall" , "Malaria_CM", "Malaria_SMA", "Malaria_Other", "Neutrophil", "Neutrophil", "Neutrophil", "Neutrophil"),
    outcome_path = c(
        "#.txt", 
        "#.txt", 
        "#.txt", 
        "#.txt",
        "#.bgen",
        "#.bgen",
        "#.bgen",
        "#.bgen"
        ),
    reverse = c("FALSE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE", "TRUE", "TRUE")
)

bim <- fread("#.bim")
names(bim) <- c("Chromosome", "MarkerName", "discard", "Position", "Allele1", "Allele2")
bim[, Allele1 := tolower(Allele1)][, Allele2 := tolower(Allele2)]
bim[, Chromosome := ifelse(Chromosome < 10, paste0("0", Chromosome), as.character(Chromosome))]

# Run the MR analysis------
myCluster <- makeCluster(4, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

mr_list <- foreach(x = 1:8) %dopar% {
  malaria_mr(
    exp_path = mr_info[[x,2]],
    exp_name = mr_info[[x,1]],
    out_path = mr_info[[x,4]],
    out_name = mr_info[[x,3]],
    reverse = mr_info[[x,5]]
  )
}

stopCluster(myCluster)

save(mr_list, file = "#.RData")

# 2S UVMR Conditional F-statistic------
fstat_fun <- function(res_df, trait_data) {
  snp_dat <- res_df[res_df$mr_keep == TRUE, ]
  snp_dat <- as.data.table(snp_dat)
  snp_dat[ , `:=` (exposure = trait_data[1], outcome = trait_data[2])]
  cols_to_keep <- c("exposure", "outcome", "SNP", "beta.exposure", "se.exposure")
  snp_dat <- snp_dat[, ..cols_to_keep]
  snp_dat$t_stat <- snp_dat$beta.exposure/snp_dat$se.exposure
  snp_dat$f_stat <- snp_dat$t_stat^2
  snp_dat[, weak_snp := ifelse(f_stat < 10, "Yes", "No")]
  return(snp_dat)
}

myCluster <- makeCluster(4, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

uvmr_fstat_list <- foreach(x = 1:48) %dopar% {
  fstat_fun(
    res_df = mr_list[[x]][[7]],
    trait_data = c(mr_list[[x]][[2]]$exposure[1], mr_list[[x]][[2]]$outcome[1])
  )
}

stopCluster(myCluster)

uvmr_fstat_all <- rbindlist(uvmr_fstat_list)
fwrite(uvmr_fstat_all, "#.txt", quote = F, row.names = F, sep = "\t", na = NA)

