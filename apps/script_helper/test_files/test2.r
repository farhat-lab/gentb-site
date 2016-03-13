##
## Testbed function for running R from Python on VM
## Constructed to be called from Rscript
## jH
##

arg <- commandArgs(trailingOnly = TRUE)

library(jsonlite)
library(randomForest)


testfunction<-function(filename){
    library(foreign)
    
    b<-read.csv(filename)
    result<-list(a=runif(1),bsum=sum(b))
    
    result<-jsonlite:::toJSON(result)
    return(result)
}

testfunction(arg)