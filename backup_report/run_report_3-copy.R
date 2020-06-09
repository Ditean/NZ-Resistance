library(knitr)
library(dplyr)
#library(ape) #needed in sweave doc
#library(ggtree) #needed in sweave doc

options(warn=-1)

makeReport<-function(reportData=NULL, drugDat = NULL, pipeDat = NULL){
  id<-reportData$BARCODE
  #copy the very original template
  con<-file.copy("./reportTemplate/resistance_script.Rnw",
                 paste0("./output/",id,"_testReport",".Rnw", collapse =""),overwrite = TRUE)
  
  #environment variables magically get passed when compiled
  Sweave2knitr(paste0("./output/",id,"_testReport",".Rnw", collapse =""),
               paste0("./output/",id,"_testReport",".Rnw"))
  con<-knit2pdf(paste0("./output/",id,"_testReport",".Rnw", collapse =""),
                compiler = 'xelatex',quiet = T)
  #moving and getting rid of extra files (just keep the pdf)
  file.remove(paste(paste0(id,"_testReport",collapse=""),c(".log",".tex",".aux"),sep=""))
  file.rename(from=paste0(id,"_testReport.pdf",collapse=""), 
              to=paste0("./output/",id,"_testReport.pdf",collapse=""))
  file.remove(paste0("./output/",id,"_testReport",".Rnw"))
}

# Testing it out with some data
#patDat<-read.csv(file="./data/patientMetaData.csv",header=T,stringsAsFactors = F)

#Note, not tidy data, which is being converted into a tiday format
#patDrugResistData<-read.csv(file = "./data/patientDrugSusceptibilityData.csv",header=T,stringsAsFactors = F)

#==============
# Patient Resistance Profile
patientResistance <- read.delim("data/resistance_profile.tsv")

# Patient Metadata
patientMetadata <- read.delim("data/patientdata.tsv")

# Pipeline Metadata
pipelineMetadata <- read.delim("data/manifest.tsv")

#Not sure how different groups will handle the cluster data, I've just opted to do something
#Simple and random that illustrates that point. I will assume the phylogeny information comes
#from some other place. Everything is generated in the sweave template as a random tree

#Generating the reports
for(i in 1:nrow(patientMetadata)){
  reportData = patientMetadata[i,]
  makeReport(reportData, filter(patientResistance,BARCODE == reportData$BARCODE), filter(pipelineMetadata,BARCODE == reportData$BARCODE))
}