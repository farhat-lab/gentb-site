# -*- coding: utf-8 -*-
# Read raw data
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# June 20, 2016

import pandas as pd

# Training Sample -- All the Mutations
data = pd.read_csv("./input/str.csv")

# Create target variable
data['y'] = (data['dr'] == "r") * 1
data.drop('dr', axis=1, inplace=True)

# List of Features to Keep in the Analysis
features = [var for var in data.columns if var != "y"]

# List subset of Features
features_small = ["SNP_CN_781687_A128G_K43R_rpsL", "SNP_N_1472359_A514C_rrs", "SNP_CN_781822_A263C_K88T_rpsL", "SNP_N_1473246_A1401G_rrs", 
				  "SNP_CN_781822_A263G_K88R_rpsL", "SNP_CN_4407809_C394A_D132Y_gid", "SNP_CN_4408156_A47C_L16R_gid", "SNP_N_1472358_C513T_rrs", "SNP_CN_4407927_T276G_E92D_gid", 
				  "SNP_N_1472751_A906G_rrs", "SNP_CN_4407934_A269C_L90R_gid", "SNP_N_1472362_C517T_rrs", "SNP_N_1472753_A908C_rrs", "SNP_CN_781822_A263T_K88M_rpsL",
				  "SNP_CN_4407832_A371G_V124A_gid", "SNP_CN_4408091_G112T_P38T_gid", "SNP_I_1473637_A21G_inter_rrs_rrl", "SNP_CN_4408094_C109T_G37R_gid", "DEL_CF_4407640_d563A_188_gid",
				  "SNP_N_1473109_T1264G_rrs", "SNP_CN_4407967_A236C_L79W_gid", "SNP_CN_4407967_A236G_L79S_gid", "SNP_CN_4407768_C435A_L145F_gid", "SNP_CN_4407995_T208G_S70R_gid",
				  "DEL_CF_4407852_d351C_117_gid", "SNP_N_1473167_T1322G_rrs", "DEL_CF_4408023_d180T_60_gid", "DEL_CF_4408116_d87G_29_gid", "SNP_CN_4408060_T143G_H48P_gid",
				  "SNP_CN_4408138_T65C_Y22C_gid", "SNP_CN_4408064_G139A_R47W_gid", "SNP_CN_4408148_C55G_A19P_gid", "SNP_CN_4407947_G256A_L86F_gid", "SNP_CN_4407916_C287A_R96L_gid",
				  "SNP_CN_4407748_A455G_L152S_gid", "SNP_N_1473343_G1498T_rrs", "SNP_N_1472337_C492T_rrs", "SNP_CN_4407985_C218G_G73A_gid", "SNP_CN_4408102_C101T_G34E_gid"]



