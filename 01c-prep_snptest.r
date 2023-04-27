# Load packages------
load_pkg <- rlang::quos(tidyverse, data.table, ggplot2, ggpubr, ggrepel, RColorBrewer, car, speedglm, foreach, doParallel, kableExtra, meta)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

# Load data------
covars <- fread("#.txt", stringsAsFactors = F)
covars <- covars[, c("genid", "neutrophil_count", "nc_log", "neutrophil_device", "assessment_center",
                     "sample_year", "sample_month", "sample_day", "sample_day_minutes", 
                     "sex", "age", paste0("PC", 1:100), "Kpop", "rs2814778_T")]
for (i in names(covars)) {
  covars[get(i) == " " | get(i) == "", (i) := NA]
}

# Create blank data.table with what snptest will use as input------
# As discussed, will run linear model in R first, then use residuals as phenotype for snptest
snptest_covars <- data.table(
    ID_1 = covars$genid, ID_2 = covars$genid, 
    missing = 0, sex = covars$sex
)

covars <- covars[!is.na(nc_log),]

# Convert categorical columns to factors------
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
covars$assessment_center <- as.factor(covars$assessment_center)

# Run a linear model for each Kpop and whole CAG, with and without accounting for Duffy allele status
lm_names <- c("res_wod_K1", "res_wod_K2", "res_wod_K3", "res_wod_K4", 
              "res_wod_K5", "res_wod_K6", "res_wod_K7", "res_wd_K1", 
              "res_wd_K2", "res_wd_K3", "res_wd_K4", "res_wd_K5", 
              "res_wd_K6", "res_wd_K7", "res_wod_all", "res_wd_all")
lm_snptest <- function(x) {
  if (grepl("wd", lm_names[x], fixed=TRUE)) {
    my_lm_vars = c("sex", "neutrophil_device", "sample_year", "sample_month",
                   "sample_day", "sample_day_minutes", "age", "assessment_center", paste0("PC", 1:100), "rs2814778_T")
  } else {
    my_lm_vars = c("sex", "neutrophil_device", "sample_year", "sample_month",
                   "sample_day", "sample_day_minutes", "age", "assessment_center", paste0("PC", 1:100))
  }
  my_lm_vars <- paste0(my_lm_vars, collapse = " + ")
  lm_formula <- formula(paste0("nc_log ~ ", my_lm_vars))
  if ( x <= 14) {
    if (x <= 7) {
      kmeans <- paste0("K", x)
    } else if (x > 7) {
      kmeans <- paste0("K", x - 7)
    }
    my_lm <- lm(lm_formula, data = covars[Kpop == kmeans,])
    lm_residuals_dt <- data.table(genid = covars[Kpop == kmeans, genid], residuals = residuals(my_lm))
    setnames(lm_residuals_dt, "residuals", lm_names[x])
  } else {
    my_lm <- lm(lm_formula, data = covars)
    lm_residuals_dt <- data.table(genid = covars[['genid']], residuals = residuals(my_lm))
    setnames(lm_residuals_dt, "residuals", lm_names[x])
  }
  return(lm_residuals_dt)
}

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

myCluster <- makeCluster(6, type = "FORK", outfile="#.log")
registerDoParallel(myCluster)

lm_list <- list()
lm_list <- foreach(x = 1:16) %dopar% {
  lm_snptest(
    x = x
  )
}
lm_list <- Reduce(function(...) merge(..., by = c('genid'), all=TRUE), lm_list)

stopCluster(myCluster)

snptest_covars <- merge(snptest_covars, lm_list, by.x = "ID_1", by.y = "genid", all = TRUE)
all_cols <- colnames(snptest_covars)
snptest_covars[, lapply(.SD, function(x) sum(!is.na(x))), .SDcols = 1:ncol(snptest_covars)]

cols <- colnames(snptest_covars)[c(5:20)]
snptest_covars[ , (cols) := lapply(.SD, signif), .SDcols = cols]

# Order to make snptest happy
all_ids <- fread("#.sample")
all_ids <- all_ids[-1,]
snptest_covars_fin <- merge(all_ids[, c("ID_1")], snptest_covars, by="ID_1", all.x=TRUE, sort=FALSE)
snptest_covars_fin[, ID_2 := ID_1]

# Final touches
snptest_covars_fin <- rbind(as.data.table(t(c("0", "0", "0", "D", rep("P", 16)))), snptest_covars_fin, use.names=FALSE)
colnames(snptest_covars_fin) <- all_cols

# Write to file
fwrite(snptest_covars_fin, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
