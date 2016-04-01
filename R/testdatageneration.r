## File to create test data for visualizer for limited mutation set tree predictions
## jH

library("jsonlite")

id <- c('inh','rif','pza', 'emb','str','eth','kan', 'cap', 'amk', 'cip', 'levo', 'oflx', 'moxi', 'gati', 'pas')
n.d<-length(id)

aa<-runif(n.d)
bb<-runif(n.d)
cc<-runif(n.d)
dd<-runif(n.d)
ee<-runif(n.d)
data<-data.frame(id,aa,bb,cc,dd,ee)
write.csv(data,file="testdata.csv",row.names=FALSE)


n.outcomes<-(2^ncol(data)) -1
mynames<-NULL
output<-list()
for(i in 0:n.outcomes){
    #cat(i,"\n")
  mynames<- c(mynames, paste("outcome",i,sep=""))
  sims<-runif(n.d)
  output[[i+1]]<- sims
}
names(output)<-mynames

result<- jsonlite:::toJSON(output)

sink("json.txt")
cat(result)
sink()



  



