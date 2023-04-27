# Load packages------
load_pkg <- rlang::quos(tidyverse, data.table, ggplot2, ggpubr, ggrepel, RColorBrewer, car)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

# Load data------
covars <- fread("#.txt", stringsAsFactors = F)
covars <- covars[, c("genid", "neutrophil_count", "nc_log", "neutrophil_device", 
                     "sample_year", "sample_month", "sample_day", "sample_day_minutes", 
                     "sex", "age", "menopause", "cob_un_regions", "Kpop", 
                     paste0("PC", 1:40),
                     "rm_astle", "BMI", "height", "weight", "home_north", 
                     "home_east", "townsendDI", "smoking_status", "drinker_status", 
                     "rs2814778_T")]
covars <- covars[!is.na(neutrophil_count), ]
for (i in names(covars)) {
  covars[get(i) == " " | get(i) == "", (i) := NA]
}
covars$menopause <- as.character(covars$menopause)
covars[is.na(menopause), menopause := ifelse(sex == 1, "male", NA)]

# Missing bar plot------
f_dowle2 = function(DT) {
  for (i in names(DT)) {
    DT[get(i) == " " | get(i) == "", (i) := NA]
    if(!i %in% c("genid", "neutrophil_count", "nc_log")) DT[, (i) := fcase(
      !is.na(get(i)) & get(i) != "prefer_not_answer" & get(i) != "-3" & get(i) != "3", 0,
      get(i) == "prefer_not_answer", 1,
      is.na(get(i)), 2,
      get(i) == "-3", 1,
      get(i) == "3", 3)]
  }
}

covars_mis <- copy(covars)
f_dowle2(covars_mis)
cols_missing <- unique(unlist(unlist(apply(covars_mis, 1, function(x) list(names(which(x == 1 | x == 2 | x == 3)))), recursive=FALSE)))
covars_mis <- covars_mis[, c("genid", "nc_log", cols_missing), with = F]

# Barplot
dt_mis_count = as.data.table(melt(covars_mis, id.vars = c("genid", "nc_log"),
                measure.vars = cols_missing) %>% 
                group_by(variable, value) %>% 
                arrange(variable, value) %>% 
                count())
dt_mis_count$variable <- as.character(dt_mis_count$variable)
dt_mis_count <- dt_mis_count[!value == 0]
dt_mis_count <- dt_mis_count[!variable %in% c("sample_day", "sample_month"),]
dt_mis_count[, pct_exp := round(n/6133*100, digits=2)]
dt_mis_count$value <- as.character(dt_mis_count$value)
dt_mis_count[value == "1", value := "Prefer not answer"]
dt_mis_count[value == "2", value := "Missing"]
dt_mis_count[value == "3", value := "Not sure"]
dt_mis_count <- dt_mis_count[order(variable, value)]
missing_barplot <- ggbarplot(dt_mis_count, "variable", "pct_exp",
   fill = "value", color = "black", xlab = "Variable name", ylab = "Value (%)", 
   label = TRUE, position = position_dodge(0.9), palette = c("firebrick", "dodgerblue", "darkgreen")) + 
   scale_x_discrete(breaks = c("BMI", "cob_un_regions", "drinker_status", "height", "menopause", "smoking_status", "townsendDI", "weight"), 
                    labels = c("Body mass index", "UN region", "Drinker status", "Height", "Menopause status", "Smoking status", "Townsend\ndepravation index", "Weight")) + 
   guides(colour = guide_legend(reverse=F)) + rotate_x_text(45)

ggsave("#.tiff", 
  missing_barplot, device="tiff", width = 210, height = 210, dpi=300, compression = "lzw", units="mm", bg="white")

# Missing model------
covars_mis[, (cols_missing) := lapply(.SD, as.factor), .SDcols = cols_missing]
missing_lm_1 <- lm(nc_log ~ menopause, data = covars_mis)
missing_lm_2 <- lm(nc_log ~ cob_un_regions, data = covars_mis)
missing_lm_3 <- lm(nc_log ~ BMI, data = covars_mis)
missing_lm_4 <- lm(nc_log ~ height, data = covars_mis)
missing_lm_5 <- lm(nc_log ~ weight, data = covars_mis)
missing_lm_6 <- lm(nc_log ~ townsendDI, data = covars_mis)
missing_lm_7 <- lm(nc_log ~ smoking_status, data = covars_mis)
missing_lm_9 <- lm(nc_log ~ drinker_status, data = covars_mis)

missing_lm_all <- lm(nc_log ~ menopause + cob_un_regions + BMI + height + weight + townsendDI + smoking_status + drinker_status, data = covars_mis)

missing_lm_list <- list(missing_lm_1, missing_lm_2, missing_lm_3, missing_lm_4, missing_lm_5, missing_lm_6, missing_lm_7, missing_lm_9, missing_lm_all)
save(missing_lm_list, file = "#.RData")

# Combine all into a data.table and curate to add as supplementary table
lm_df <- data.table(summary(missing_lm_list[[1]])[[4]], keep.rownames = TRUE)
for (x in 2:8) {
  lm_df <- rbind(lm_df, data.table(summary(missing_lm_list[[x]])[[4]], keep.rownames = TRUE))
}

lm_df <- lm_df[rn != "(Intercept)", ]
lm_df[, type := fcase(
  grepl("1", rn), "Prefer not answer", 
  grepl("2", rn), "Missing", 
  grepl("3", rn), "Not sure"
)]
lm_df$rn <- substr(lm_df$rn, 1, nchar(lm_df$rn) - 1)
lm_df <- merge(lm_df, dt_mis_count, by.x = c("rn", "type"), by.y = c("variable", "value"))
lm_df$outcome <- "log_neutrophil_count"
lm_df <- lm_df[, c("rn", "outcome", "type", "Estimate", "Std. Error", "Pr(>|t|)", "n")]
setnames(lm_df, c("rn", "outcome", "type", "Estimate", "Std. Error", "Pr(>|t|)", "n"), 
  c("exposure", "outcome", "type", "beta", "se", "pval", "n"))
fwrite(lm_df, "#.txt", quote=F, row.names=F, na=NA, sep="\t")

