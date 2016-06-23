# -*- coding: utf-8 -*-
# Read raw data
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# May 15, 2016

setwd("t://share//jimmy//gentb-site//R//Neural_Network")
print(getwd())

load("../eth_finalpredict.RData")
write.csv(selected, "./input/short_eth.csv")
write.csv(drugg.full, "./input/eth.csv", row.names=FALSE)

load("../str_finalpredict.RData")
write.csv(selected, "./input/short_str.csv")
write.csv(drugg.full, "./input/str.csv", row.names=FALSE)

load("../emb_finalpredict.RData")
write.csv(selected, "./input/short_emb.csv")
write.csv(drugg.full, "./input/emb.csv", row.names=FALSE)
