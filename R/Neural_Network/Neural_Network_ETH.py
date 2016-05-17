# -*- coding: utf-8 -*-
# 2-hidden Layer Neural Network
# Compute ROC, Sensitivity, Marginal Effects
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# May 16, 2016

# Uncomment to Run on GPU
#from sknn.platform import gpu32
import os
import numpy as np

# Change to Current Directory
os.chdir(r'\gentb-site\R\Neural_Network')

# Load Data
execfile('Load_Data.py')

# Helper Functions
execfile('Helper.py')

## Collect Feature Names
predictors = np.chararray((len(features), 1), itemsize=50)
for i in range(len(features)):
    predictors[i,0] = features[i]

## Validation Set
X = data[features]
y = data["y"] 

## Output Containers and Boostrap Paramters
n_boot = 100
gof_measures = np.zeros((3, (n_boot+2)), dtype=np.float32)
marg_effects = np.zeros((len(features), (n_boot+3)), dtype=np.float32)

## Grid Search for Meta-Parameters
Best = meta(X, y)

## Bootstrap Marginal Effects w/ Standard Errors
for i in range(n_boot):
    boot(i, Best, X, y)

## Compute Mean and Standard Deviation of Bootstraped Marginal Effects and Measures of Fit
marg_effects[:, n_boot] = np.mean(marg_effects[:,range(n_boot)], axis=1) 
marg_effects[:, (n_boot+1)] = np.std(marg_effects[:,range(n_boot)], axis=1) 
marg_effects[:, (n_boot+2)] = marg_effects[:, n_boot] / marg_effects[:, (n_boot+1)]

gof_measures[:, n_boot] = np.mean(gof_measures[:,range(n_boot)], axis=1)
gof_measures[:, (n_boot+1)] = np.std(gof_measures[:,range(n_boot)], axis=1)

## Export Results
toout = np.concatenate([predictors, marg_effects], axis=1)
ind = np.abs(np.double(toout[:, (n_boot+3)])).argsort()[::-1]
sortedout = toout[ind]

np.savetxt("./marginal_effects.csv", sortedout, fmt="%s", delimiter=',')
np.savetxt("./GofF.csv", gof_measures, fmt="%s", delimiter=',')
