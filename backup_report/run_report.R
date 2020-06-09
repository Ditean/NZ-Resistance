library(knitr)
library(dplyr)

options(warn=-1)

makeReport<-function(pInfo = NULL,drugDat = NULL){
  id<-pInfo$ID
  
  con<-file.copy("/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/reportTemplate/Resistance_script.Rnw",
                 paste0("/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/output/",id,"_testReport",".Rnw", collapse = ""), overwrite = TRUE)
  
  Sweave2knitr(paste0("/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/output/",id,"_testReport",".Rnw", collapse = ""),
               paste0("/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/output/",id,"_testReport",".Rnw"))
  con<-knit2pdf(paste0("/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/ouput/",id,"_testReport",".Rnw", collapse = ""),
                compiler = 'xelatex', quiet = TRUE)
  
  file.remove(paste(paste0(id,"_testReport", collapse=""),c(".log",".tex",".aux"),sep=""))
  file.rename(from=paste0(id,"_testReport.pdf",collapse=""),
              to=paste0("/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/output/",id,"_testReport.pdf",collapse=""))
  file.remove(paste0("/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/output/",id,"_testReport",".Rnw"))
}

# Feed in data
patDat<-read.csv(file="/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/data/patientMetaData.csv",header = TRUE,stringsAsFactors = FALSE)

patDrugResistData<-read.csv(file = "/Users/tayjo75p/Documents/manuscripts/resistance_pipeline/report/data/patientDrugSusceptibilityData.csv",header=TRUE,stringsAsFactors = FALSE)

# Generate report
for(i in 1:nrow(patDat)){
  pInfo = patDat[i,]
  makeReport(pInfo,filter(patDrugResistData, ID == pInfo$ID))
}