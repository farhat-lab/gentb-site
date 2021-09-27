#!/usr/bin/env Rscript

## Minimal prediction routine with all possible preprocessing
## Configured for calling by Rscript
## Call as:
## Rscript TBpredict.R "inputfile.csv" "inputfile2.csv"

#import arguments passed to Rscript
arg <- commandArgs(trailingOnly = TRUE)

#add options (for future script version)
options <- commandArgs(trailingOnly = FALSE)
prefix <- "--file="

#source script name
script.name <- sub(prefix, "", options[grep(prefix, options)])
script.basename <- dirname(script.name)

#source Rlibray paths relative to this script's location
libs <- file.path(script.basename, '../data/Rlib')
data_dir <- file.path(script.basename, '../data/predict_rdata')

#import libraries needed for prediction
library(foreign)
library(jsonlite, quietly=TRUE, lib.loc=libs)
suppressPackageStartupMessages(library("randomForest", lib.loc=libs))

#start defining the predict function 
predictfunction<-function(filename1, filename2){
  set.seed(5414)
  
##################################################################
########### First predict using RF1.0 to 12 drugs ################
##################################################################

  druglist <- c('inh','rif','emb','str','eth','kan', 'cap', 'amk', 'cip', 'levo', 'oflx', 'pas') 
  greplist <- c('inhA|katG|embB|ahpC|ini|kasA|mabA|ndh|oxyR', 'rpoB', 'emb|ini', 'rpsL|gid|rrs', 'ethA|inhA|fabG1', 'rrs|tlyA', 'tlyA|rrs|rrl', 'tlyA|rrs|rrl', 'gyr', 'gyr', 'gyr','thyA')
  
  input <- read.csv(filename1, header=TRUE)
  result<-matrix(NA,nrow=1,ncol=5, dimnames=list(c(), c('strain','drug','probability of resistance','False negative percent', 'False positive percent')))
  important<-vector("list", nrow(input))
  names(important)<-input[,1]
  other<-vector("list",nrow(input))
  names(other)<-input[,1]
  distmatrix<-matrix(NA,nrow=1,ncol=12,dimnames=list(c(), druglist))
  
  for ( j in 1:nrow(input)){
    strain <- input[j,]
    resultperstrain <- matrix(NA, nrow=length(druglist), ncol=5)
    important_strain<-matrix(NA,nrow=5,ncol=length(druglist),dimnames=list(c(), druglist))
    other_strain<-matrix(NA, nrow=5,ncol=length(druglist),dimnames=list(c(),druglist))
    fordiststrain<-c()
    
    for (i in 1:length(druglist)){
      drug <- druglist[i]
      
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


  result<-result[-1,] #removing empty first line
  distmatrix<-distmatrix[-1,] #removing empty first line
  if (is.vector(distmatrix)) {
    order=1
  } else {
    order<-hclust(dist(distmatrix))$order
    result2<-rep(NA,ncol(result))
    for (j in order) {
      result2<-rbind(result2, result[((j-1)*13+1):((j-1)*13+13),])

    }
    result<-result2[-1,] #removing empty first line and replacing unordered matrix
  }
  
  
##################################################################
###### Second predict using RF2.0 to Pyrazinamide  ###############
##################################################################

  # import pyrazinamide prediction matrix
  strain <- read.csv(filename2, header=TRUE)

  #prepare output for variants
  important_strain<-matrix(NA,nrow=5,ncol=1,dimnames=list(c(), 'pza'))
  important_PZA<-vector("list", nrow(strain))
  names(important_PZA)<-strain[,1]

  #prepare output for RandomForest probability
  result_PZA<-matrix(NA,nrow=1,ncol=5, dimnames=list(c(), c('strain','drug','probability of resistance','False negative percent', 'False positive percent')))

  #load RandomForest object
  load("data/predict_rdata/pza_finalpredict_v2_0.RData")

  #perform prediction and write strain ID, drug, and probability to result
  Valid <- predict(drugg.full.rf, strain, type='prob', norm.votes=TRUE, predict.all=FALSE)
  result_PZA[1,] <- c(as.vector(strain[1,1]), 'pza', round(Valid[1,1],3), '16.1', '7.6')

  #write the variant used for prediction to object
  imp<-colnames(strain)[which(strain[1,]==1)]
  
  #in case no variants were detected write NA's 
  if(length(imp) == 0) {imp <- c(NA,NA,NA,NA,NA)}
  important_strain[1:length(imp),1]<-imp
  important_PZA[[1]]<-important_strain

  #bind the probabilities from RF1.0 and RF2.0
  result_both_predict <- rbind(result, result_PZA)

  #bind important variant lists from RF1.0 and RF2.0
  important_both_predict <- list(mapply(c, important, important_PZA, SIMPLIFY=FALSE))

  
  #make a list of RF1.0 and RF2.0 probabilities and variants
  listOfpredictions <- list(result_both_predict, important_both_predict, other)

  # Save JSON file
  file_noext <- substr(filename2, 1, nchar(filename2) - 24)
  cat(toJSON(listOfpredictions, pretty = TRUE), "\n", file = paste0(file_noext, ".matrix.json"))
  
}

#Call prediction function on both input arguments
suppressWarnings(predictfunction(arg[1], arg[2]))
