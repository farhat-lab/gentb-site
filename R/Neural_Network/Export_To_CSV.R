
setwd("./Neural_Network")
load("../R/eth_finalpredict.RData")
getwd()

drugg.full["y"] = 0
drugg.full$y = drugg.full$dr == "r"

write.csv(drugg.full, "./eth.csv")