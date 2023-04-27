library(data.table)

bim <- fread("#")
names(bim) <- c("Chromosome", "MarkerName", "discard", "Position", "Allele1", "Allele2")
bim[, Allele1 := tolower(Allele1)][, Allele2 := tolower(Allele2)]
bim[, Chromosome := ifelse(Chromosome < 10, paste0("0", Chromosome), as.character(Chromosome))]

# bolt p<5e-5
bolt <- fread("#.bgen")
fwrite(bolt[P_BOLT_LMM_INF < 5e-5,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Load malariagen overall
overall_sm <- fread("#.txt")
overall_sm_merged_A1 <- merge(overall_sm, bim, by=c("Chromosome", "Position", "Allele1"))
overall_sm_merged_A1 <- overall_sm_merged_A1[, c("Chromosome", "Position", "MarkerName.y", "Allele1", "Allele2.x", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction")]
setnames(overall_sm_merged_A1, c("MarkerName.y", "Allele2.x"), c("MarkerName", "Allele2"))

overall_sm_merged_A2 <- merge(overall_sm, bim, by=c("Chromosome", "Position", "Allele2"))
overall_sm_merged_A2 <- overall_sm_merged_A2[, c("Chromosome", "Position", "MarkerName.y", "Allele1.x", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction")]
setnames(overall_sm_merged_A2, c("MarkerName.y", "Allele1.x"), c("MarkerName", "Allele1"))

overall_sm_merged <- rbind(overall_sm_merged_A1, overall_sm_merged_A2)
fwrite(overall_sm_merged, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(overall_sm_merged[`P-value` < 5e-5,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(overall_sm_merged[`P-value` < 5e-8,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Load malariagen CM
cm_sm <- fread("#.txt")
cm_sm_merged_A1 <- merge(cm_sm, bim, by=c("Chromosome", "Position", "Allele1"))
cm_sm_merged_A1 <- cm_sm_merged_A1[, c("Chromosome", "Position", "MarkerName.y", "Allele1", "Allele2.x", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction")]
setnames(cm_sm_merged_A1, c("MarkerName.y", "Allele2.x"), c("MarkerName", "Allele2"))

cm_sm_merged_A2 <- merge(cm_sm, bim, by=c("Chromosome", "Position", "Allele2"))
cm_sm_merged_A2 <- cm_sm_merged_A2[, c("Chromosome", "Position", "MarkerName.y", "Allele1.x", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction")]
setnames(cm_sm_merged_A2, c("MarkerName.y", "Allele1.x"), c("MarkerName", "Allele1"))

cm_sm_merged <- rbind(cm_sm_merged_A1, cm_sm_merged_A2)
fwrite(cm_sm_merged, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(cm_sm_merged[`P-value` < 5e-5,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(cm_sm_merged[`P-value` < 5e-8,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Load malariagen SMA
sma_sm <- fread("#.txt")
sma_sm_merged_A1 <- merge(sma_sm, bim, by=c("Chromosome", "Position", "Allele1"))
sma_sm_merged_A1 <- sma_sm_merged_A1[, c("Chromosome", "Position", "MarkerName.y", "Allele1", "Allele2.x", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction")]
setnames(sma_sm_merged_A1, c("MarkerName.y", "Allele2.x"), c("MarkerName", "Allele2"))

sma_sm_merged_A2 <- merge(sma_sm, bim, by=c("Chromosome", "Position", "Allele2"))
sma_sm_merged_A2 <- sma_sm_merged_A2[, c("Chromosome", "Position", "MarkerName.y", "Allele1.x", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction")]
setnames(sma_sm_merged_A2, c("MarkerName.y", "Allele1.x"), c("MarkerName", "Allele1"))

sma_sm_merged <- rbind(sma_sm_merged_A1, sma_sm_merged_A2)
fwrite(sma_sm_merged, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(sma_sm_merged[`P-value` < 5e-5,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(sma_sm_merged[`P-value` < 5e-8,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")

# Load malariagen OTHER
other_sm <- fread("#.txt")
other_sm_merged_A1 <- merge(other_sm, bim, by=c("Chromosome", "Position", "Allele1"))
other_sm_merged_A1 <- other_sm_merged_A1[, c("Chromosome", "Position", "MarkerName.y", "Allele1", "Allele2.x", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction")]
setnames(other_sm_merged_A1, c("MarkerName.y", "Allele2.x"), c("MarkerName", "Allele2"))

other_sm_merged_A2 <- merge(other_sm, bim, by=c("Chromosome", "Position", "Allele2"))
other_sm_merged_A2 <- other_sm_merged_A2[, c("Chromosome", "Position", "MarkerName.y", "Allele1.x", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction")]
setnames(other_sm_merged_A2, c("MarkerName.y", "Allele1.x"), c("MarkerName", "Allele1"))

other_sm_merged <- rbind(other_sm_merged_A1, other_sm_merged_A2)
fwrite(other_sm_merged, "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(other_sm_merged[`P-value` < 5e-5,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")
fwrite(other_sm_merged[`P-value` < 5e-8,], "#.txt", quote=F, row.names=F, na=NA, sep="\t")
