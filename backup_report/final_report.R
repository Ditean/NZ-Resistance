library(knitr)
library(dplyr)

options(warn = -1)

makeReport <- function(repData = NULL, repDrug = NULL, repMeta = NULL, repSNP = NULL){
  barcode <- repData$BARCODE

  repTmp <- file.copy("./reportTemplate/tutorial.Rnw",
                      paste0("./output/",barcode,"_report",".Rnw", collapse = ""), overwrite = TRUE)

  Sweave2knitr(paste0("./output/",barcode,"_report",".Rnw", collapse = ""),
               paste0("./output/",barcode,"_report",".Rnw"))
  repTmp <- knit2pdf(paste0("./output/",barcode,"_report",".Rnw", collapse = ""),
                     compiler = 'xelatex', quiet = TRUE)

  file.remove(paste(paste0(barcode,"_report", collapse = ""), c(".log", ".tex", ".aux"), sep = ""))
  file.rename(from=paste0(barcode,"_report.pdf", collapse = ""),
              to=paste0("./output/", barcode, "_report.pdf", collapse=""))
  file.remove(paste0("./output/", barcode, "_report", ".Rnw"))
}

# Detected resistance file
crudeDrug <- read.delim("./data/resistance_profile.tsv", stringsAsFactors = FALSE)

# Patient Metadata
crudeMeta <- read.delim("./data/patientdata.tsv", stringsAsFactors = FALSE)

# Pipeline information
crudePipe <- read.delim("./data/manifest.tsv", stringsAsFactors = FALSE)

crudeSNP <- read.delim("./data/genomic_profile.tsv", stringsAsFactors = FALSE)

# Report automation
for(i in 1:nrow(crudeMeta)){
  repData = crudeMeta[i,]
  makeReport(repData, filter(crudeDrug,BARCODE == repData$BARCODE), filter(crudePipe,BARCODE == repData$BARCODE), filter(crudeSNP,BARCODE == repData$BARCODE))
}
