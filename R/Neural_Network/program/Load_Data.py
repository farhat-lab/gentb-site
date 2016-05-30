# -*- coding: utf-8 -*-
# Read raw data
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# May 15, 2016

import pandas as pd

# Training Sample -- All the Mutations
data = pd.read_csv("./input/eth.csv")

# Create target variable
data['y'] = (data['dr'] == "r") * 1
data.drop('dr', axis=1, inplace=True)

# List of Features to Keep in the Analysis
features = [var for var in data.columns if var != "y"]

# List subset of Features
features_small = ["SNP_P_1673425_C15T_promoter_fabG1.inhA", "SNP_CN_4326333_C1141G_A381P_ethA", "SNP_CN_4326116_G1358A_T453I_ethA", "SNP_CN_1674481_T280G_S94A_inhA",
                  "SNP_CZ_4326714_G760A_Q254._ethA", "SNP_CN_1674263_T62C_I21T_inhA", "SNP_CN_4327416_C58A_A20S_ethA", "DEL_CF_4326184_d1290G_430_ethA", 
                  "SNP_CN_4327380_A94C_Y32D_ethA", "SNP_CN_1674434_T233G_V78G_inhA", "INS_CF_4326141_i1333C_445_ethA", "SNP_CZ_4326600_G874A_R292._ethA",
                  "SNP_CN_4326713_T761G_Q254P_ethA", "SNP_CN_4326305_G1169A_S390F_ethA", "SNP_P_1673423_G17T_promoter_fabG1.inhA", "INS_CF_4326722_i752C_251_ethA", 
                  "SNP_CN_1673449_A10C_T4P_fabG1", "SNP_CN_4327311_A163G_S55P_ethA", "SNP_CZ_4326278_G1196T_S399._ethA", "SNP_CZ_4327148_C326T_W109._ethA"]


