# -*- coding: utf-8 -*-
# Read raw data
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# June 20, 2016

import pandas as pd

# Training Sample -- All the Mutations
data = pd.read_csv("./input/emb.csv")

# Create target variable
data['y'] = (data['dr'] == "r") * 1
data.drop('dr', axis=1, inplace=True)

# List of Features to Keep in the Analysis
features = [var for var in data.columns if var != "y"]

# List subset of Features
features_small = ["SNP_CN_4247429_A916G_M306V_embB", "SNP_CN_4247431_G918A_M306I_embB", "SNP_CN_4247431_G918C_M306I_embB", "SNP_CN_4247730_G1217C_G406A_embB", 
				  "SNP_CN_4248003_A1490G_Q497R_embB", "SNP_CN_4249518_A3005G_H1002R_embB", "SNP_CN_409569_G208A_A70T_iniB", "SNP_CN_4247729_G1216A_G406S_embB", 
				  "SNP_CN_4247431_G918T_M306I_embB", "SNP_CN_4247429_A916C_M306L_embB", "SNP_P_4243222_C11A_promoter_embA.embB", "SNP_CN_4247574_A1061C_D354A_embB", 
				  "SNP_CN_4247495_G982T_D328Y_embB", "SNP_CN_4249583_G3070A_D1024N_embB", "SNP_CN_4243392_A160G_N54D_embA", "SNP_P_4243225_C8T_promoter_embA.embB", 
				  "SNP_CN_4242182_G2320T_A774S_embC", "SNP_CN_4247729_G1216T_G406C_embB"]