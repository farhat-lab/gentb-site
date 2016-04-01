##
## Second testbed function for running R from Python on VM
## Constructed to be called from Rscript
## jH
##

arg <- commandArgs(trailingOnly = TRUE)

library(jsonlite)
library(randomForest)
library(foreign)

testfunction<-function(filename){

  prddata <- read.csv(arg)
  load("outRF.Rdata")
  prd<- predict(object=out, newdata=prddata, type="prob", norm.votes=TRUE, predict.all=FALSE)
  keep<-list(r=prd[1],s=prd[2])
  result<- jsonlite:::toJSON(keep)     # rjson does not format json correctly for a list of lists
  return(result)
}

testfunction(arg)

