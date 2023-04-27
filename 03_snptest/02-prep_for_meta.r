# Load relevant libraries------
list_of_packages <- c(
  "data.table", "tidyverse",
  "ggplot2", "ggpubr", "ggrepel", "RColorBrewer",
  "ggforestplot"
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

# Filter each chromosome for each kpop
# Replace SNP beta, se, and pvalue with NA if info < 0.3 or/and MAC < 10

my_snptest_dirs <- list.files("#", pattern="K.*", full.names=T)
my_chrs <- paste0(sprintf("chr%02d", 1:22),".txt")


for (x in my_snptest_dirs) {
    for(y in my_chrs) {
        print(paste0("Now at dir ", x, " and file ", y))
        myfile <- paste0(x, "/", y)
        mydat <- fread(myfile)
        mydat[, mac := all_maf * all_total]
        mydat <- mydat[mac > 10 & info > 0.3,]
        fwrite(mydat, paste0(x, "/", y, ".filtered"), quote=F, row.names=F, na=NA, sep="\t")
    }
}

