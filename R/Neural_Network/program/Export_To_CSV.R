# -*- coding: utf-8 -*-
# Read raw data
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# May 15, 2016

setwd("../")
print(getwd())
load("../eth_finalpredict.RData")

write.csv(drugg.full, "./input/eth.csv", row.names=FALSE)
