## Minimal prediction routine with all possible preprocessing
## Configured for calling by Rscript
## Call as:
## Rscript TBpredict.R "testinputfile.csv"


arg <- commandArgs(trailingOnly = TRUE)

aa <- Sys.time()

local({r<-getOption("repos")
	 r["CRAN"]<-"http://cran.fhcrc.org/"
	 options(repos=r)})

#install.packages("jsonlite", quiet=TRUE)
library("foreign")
suppressPackageStartupMessages(library(jsonlite, quietly=TRUE, lib="/www/gentb.hms.harvard.edu/support/R.2.15.3.lib"))
suppressPackageStartupMessages(library(randomForest, quietly=TRUE, lib="/www/gentb.hms.harvard.edu/support/R.2.15.3.lib"))


predictfunction<-function(filename){
  
  set.seed(5414)
  
  druglist <- c('inh','rif','pza', 'emb','str','eth','kan', 'cap', 'amk', 'cip', 'levo', 'oflx', 'pas') 
  greplist <- c('inhA|katG|embB|ahpC|ini|kasA|mabA|ndh|oxyR', 'rpoB', 'pncA', 'emb|ini', 'rpsL|gid|rrs', 'ethA|inhA|fabG1', 'rrs|tlyA', 'tlyA|rrs|rrl', 'tlyA|rrs|rrl', 'gyr', 'gyr', 'gyr','thyA')
  
  input <- read.csv(filename, header=TRUE)
  result<-matrix(NA,nrow=1,ncol=5, dimnames=list(c(), c('strain','drug','probability of resistance','False negative percent', 'False positive percent')))
  important<-vector("list", nrow(input))
  names(important)<-input[,1]
  other<-vector("list",nrow(input))
  names(other)<-input[,1]
  
  for ( j in 1:nrow(input)){
    strain <- input[j,]
    resultperstrain <- matrix(NA, nrow=length(druglist), ncol=5)
    important_strain<-matrix(NA,nrow=5,ncol=length(druglist),dimnames=list(c(), druglist))
    other_strain<-matrix(NA, nrow=5,ncol=length(druglist),dimnames=list(c(),druglist))
    
    for (i in 1:length(druglist)){
      drug <- druglist[i]
      load(paste("/groups/murray/run_pipeline/bin/", drug, "_finalpredict.RData", sep="")) #contains for each drug:
      #drugg.full   :   matrix of all data
      #drugg.full.rf:   a RF classifier 
      #eR           :   OOB false negative rate (For resistance classification)
      #eS           :   OOB false positive rate (For resistance classification)
      #selected     :   important variables for resistance prediction in order of their prediction importance
      
      #filtering nonsyn mutations in genes of interest for each drug, it's not really necessary but may add speed
      muts <- grep(greplist[i],colnames(strain))
      nonsynmuts <- grep('CS', invert=TRUE, colnames(strain))
      vars <- strain[,intersect(nonsynmuts, muts)]

      Valid <- predict(drugg.full.rf, vars, type="prob", norm.votes=TRUE, predict.all=FALSE)
      resultperstrain[i,] <- c(as.vector(strain[1,1]), drug, round(Valid[1,1],3), round(eR,2), round(eS,2))
      
      imp<-intersect(selected, colnames(strain)[which(strain[1,]==1)])[1:5]
      important_strain[1:length(imp),i]<-imp
      oth<-setdiff(intersect(selected, colnames(strain)[which(strain[1,]==1)]), colnames(strain)[which(strain[1,]==1)])[1:5]
      other_strain[1:length(oth),i]<-oth
      
    }
    result<-rbind(result, resultperstrain)
    important[[j]]<-important_strain
    other[[j]]<-other_strain
  }
  #print(length(warnings()))
  #warnings()
  #bb <- Sys.time()
  #print("Time to completion:")
  #print(bb-aa)
  result<-result[-1,] #removing empty first line
  displayresult<- jsonlite:::toJSON(list(result, important, other))
  return(displayresult)
  #return(list(result, important, other))
}

suppressWarnings(predictfunction(arg))
