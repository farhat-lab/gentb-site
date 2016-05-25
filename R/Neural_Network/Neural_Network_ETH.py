# -*- coding: utf-8 -*-
# 2-hidden Layer Neural Network
# Compute ROC, Sensitivity, Marginal Effects
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# May 25, 2016

# Uncomment to Run on GPU
#from sknn.platform import gpu32
import os
import numpy as np

# Change to Current Directory
os.chdir(r'\\mon-jroyer2\users\jroyer\desktop\dna\gentb-site\R\Neural_Network')

# Load Data
execfile('Load_Data.py')

# Helper Functions
execfile('Helper.py')

## Collect Feature Names

## Trigger Small Set of Features
# features = features_small

predictors = np.chararray((len(features), 1), itemsize=50)
for i in range(len(features)):
    predictors[i,0] = features[i]

## Validation Set
X = data[features]
y = data["y"] 

## Output Containers and Boostrap Paramters
n_boot = 250
gof_measures_NN = np.zeros((4, (n_boot+2)), dtype=np.float32)
marg_effects_NN = np.zeros((len(features), (n_boot+3)), dtype=np.float32)
gof_measures_rf = np.zeros((4, (n_boot+2)), dtype=np.float32)
marg_effects_rf = np.zeros((len(features), (n_boot+3)), dtype=np.float32)

## Grid Search for Meta-Parameters
## Need to Define Activation for Neural - Network.
## Supervised: Rectifier or ExpLin
Classifiers = meta(X, y, u'ExpLin')
NN = Classifiers[0]
rf = Classifiers[1]

## Extract Number of Units for AutoEncoder
units = []
for j in NN.layers:
    units.append(j.units)

## Run AutoEncoder
#myweights = auto(X, u'ExpLin', units)
#
### Explore Predicted X
#new_X = myweights.predict_proba(X)
#x_tild = []
#for i in new_X:
#    x_tild.append(i[:,1])
#x_tild = np.array(x_tild).T
#
### Transfer Initial Weights
#wgt = myweights.get_parameters()
wgt = NN.get_parameters()
#wgt[1] = wgt1[1]

## Bootstrap Marginal Effects w/ Standard Errors
for i in range(n_boot):
    boot(i, NN, X, y, marg_effects_NN, gof_measures_NN, "Neural Network", wgt)
    boot(i, rf, X, y, marg_effects_rf, gof_measures_rf, "Random Forest", wgt)

out_NN = margeffects(marg_effects_NN, gof_measures_NN, n_boot, predictors)
out_rf = margeffects(marg_effects_rf, gof_measures_rf, n_boot, predictors)

## Output Variable Importance Random Forest
importances = rf.feature_importances_
std = np.std([tree.feature_importances_ for tree in rf.estimators_],
             axis=0)
indices = np.argsort(importances)[::-1]
import_features = np.concatenate((predictors[indices], importances[indices].reshape(len(features),1)), axis=1)

np.savetxt("./marginal_effects_NN.csv", out_NN, fmt="%s", delimiter=',')
np.savetxt("./GofF_NN.csv", gof_measures_NN, fmt="%s", delimiter=',')
np.savetxt("./marginal_effects_rf.csv", out_rf, fmt="%s", delimiter=',')
np.savetxt("./GofF_rf.csv", gof_measures_rf, fmt="%s", delimiter=',')
np.savetxt("./importance_rf.csv", import_features, fmt="%s", delimiter=',')

