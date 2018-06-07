## Minimal prediction routine with all possible preprocessing
## Configured for calling by Rscript
## Call as:
## Rscript TBpredict.R "testinputfile.csv"

arg <- commandArgs(trailingOnly = TRUE)

options <- commandArgs(trailingOnly = FALSE)
prefix <- "--file="
script.name <- sub(prefix, "", options[grep(prefix, options)])
script.basename <- dirname(script.name)

# Set the location of the libs to relative to this script's location
libs <- file.path(script.basename, '../data/Rlib')

data_dir <- file.path(script.basename, '../data/predict_rdata')
#aa <- Sys.time()

library(foreign)
library(jsonlite, quietly=TRUE, lib.loc=libs)
suppressPackageStartupMessages(library("randomForest", lib.loc=libs))

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
  distmatrix<-matrix(NA,nrow=1,ncol=13,dimnames=list(c(), druglist))
  
  for ( j in 1:nrow(input)){
    strain <- input[j,]
    resultperstrain <- matrix(NA, nrow=length(druglist), ncol=5)
    important_strain<-matrix(NA,nrow=5,ncol=length(druglist),dimnames=list(c(), druglist))
    other_strain<-matrix(NA, nrow=5,ncol=length(druglist),dimnames=list(c(),druglist))
    fordiststrain<-c()
    
    for (i in 1:length(druglist)){
      drug <- druglist[i]
      #print(drug)
      load(paste(data_dir, "/", drug, "_finalpredict.RData", sep="")) #contains for each drug:
      # Can't get what's in RF
      #cat(toJSON(drugg.full.rf, pretty = TRUE), "\n", file = paste0(drug, "-drugg.full.rf.json"))

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
      fordiststrain<-c(fordiststrain, round(Valid[1,1],3))
      
      imp<-intersect(selected, colnames(vars)[which(vars[1,]==1)])[1:5]
      important_strain[1:length(imp),i]<-imp
      oth<-setdiff(colnames(vars)[which(vars[1,]==1)], intersect(selected, colnames(vars)[which(vars[1,]==1)]))[1:5]
      other_strain[1:length(oth),i]<-oth
      
    }
    result<-rbind(result, resultperstrain)
    important[[j]]<-important_strain
    other[[j]]<-other_strain
    distmatrix<-rbind(distmatrix,fordiststrain)
  }
  #print(length(warnings()))
  #warnings()
  #bb <- Sys.time()
  #print("Time to completion:")
  #print(bb-aa)
  result<-result[-1,] #removing empty first line
  distmatrix<-distmatrix[-1,] #removing empty first line
  if (is.vector(distmatrix)) {
    order=1
  } else {
    order<-hclust(dist(distmatrix))$order
    result2<-rep(NA,ncol(result))
    for (j in order) {
      result2<-rbind(result2, result[((j-1)*13+1):((j-1)*13+13),])
      #print((j-1)*13+1)
      #print((j-1)*13+13)
    }
    result<-result2[-1,] #removing empty first line and replacing unordered matrix
  }
  #displayresult<- jsonlite:::toJSON(result)
  #return(displayresult)
  l <- list(result, important, other)
  ## Save JSON file
  file_noext <- substr(filename, 1, nchar(filename) - 4)
  cat(toJSON(l, pretty = TRUE), "\n", file = paste0(file_noext, ".json"))
  
  #return(h)
  #return(important)
  #return(other)
}

suppressWarnings(predictfunction(arg))


# l <- predictfunction("matrix.csv")
# h <- create_heatmap("matrix.json")
# h
# 

