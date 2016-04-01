## Minimal prediction routine with all possible preprocessing
## Configured for calling by Rscript
## Call as:
## Rscript TBpredict.R "testinputfile.csv"


arg <- commandArgs(trailingOnly = TRUE)

aa <- Sys.time()

library(foreign)
library(jsonlite, quietly=TRUE)
suppressPackageStartupMessages(library("randomForest"))
library(htmlwidgets, quietly=TRUE)
library(reshape2, quietly=TRUE)
library(d3heatmap, quietly=TRUE)
library(colorRamps, quietly=TRUE)

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
      load(paste(drug, "_finalpredict.RData", sep="")) #contains for each drug:
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
  #displayresult<- jsonlite:::toJSON(result)
  #return(displayresult)
  l <- list(result, important, other)
  ## Save JSON file
  file_noext <- substr(basename(filename), 1, nchar(filename) - 4)
  cat(toJSON(l, pretty = TRUE), "\n", file = paste0(file_noext, ".json"))
  return(l)
  #return(important)
  #return(other)
}

suppressWarnings(predictfunction(arg))

create_heatmap <- function(file_json = "matrix.json") {
  l <- fromJSON(txt = readLines(file_json))
  d <- data.frame(l[[1]])[, 1:3]
  names(d) <- c("strain", "drug","probability of resistance")
  dd <- dcast(d, drug ~ strain)
  rownames(dd) <- dd$drug
  for (j in 2:ncol(dd))
    dd[, j] <- as.numeric(dd[, j])
  dd <- dd[, -1]
  dd <- t(dd)
  dd <- as.data.frame(dd)
  druglist <- c('inh','rif','pza', 'emb','str','eth','kan', 'cap', 'amk', 'cip', 'levo', 'oflx', 'pas') 
  dd <- dd[, druglist]
  return(d3heatmap(dd, dendrogram = NULL, Rowv = FALSE, Colv = FALSE, colors = colorRamps::matlab.like(10)))
}

# l <- predictfunction("matrix.csv")
# h <- create_heatmap("matrix.json")
# h
# htmlwidgets::saveWidget(h, "matrix_heatmap.html")
