# -*- coding: utf-8 -*-
# 2-hidden Layer Neural Network
# Compute ROC, Sensitivity, Marginal Effects
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# May 29, 2016

# Uncomment to Run on GPU
#from sknn.platform import gpu32
import os
import numpy as np
import pickle


# Change to Current Directory
os.chdir(r'T:/Share/Jimmy/gentb-site/R/Neural_Network')

# Load Data
execfile('./program/Load_Data_ETH.py')

# Helper Functions
execfile('./program/Helper.py')

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
#NN = find_meta_parameters(X, y, "NN", act=u'ExpLin')
#rf = find_meta_parameters(X, y, "RF")

## Save Trained Classifier
#pickle.dump(NN, open('./output/Classifiers/NN_ETH.pkl', 'wb'))
#pickle.dump(rf, open('./output/Classifiers/rf_ETH.pkl', 'wb'))
#wgtandbias = NN.get_parameters()
#np.savez('./output/Classifiers/wgtandbias', wgtandbias)

## Load Trained Classifier (Include weights and biases)
NN = pickle.load(open('./output/Classifiers/NN_ETH.pkl', 'rb'))
rf = pickle.load(open('./output/Classifiers/rf_ETH.pkl', 'rb'))

#with np.load('./output/Classifiers/wgtandbias.npz') as f:
#     param_values = [f['arr_%d' % i] for i in range(len(f.files))]
#importwgt = param_values[0]
#wgtandbias = []
#lyr = collections.namedtuple('Parameters', 'weights biases layer')
#it = 0
#for i in NN.layers:
#    lname = i.name
#    wgt = importwgt[it][0]
#    bias = importwgt[it][1]
#    wgtandbias.append(lyr(layer=lname, weights=wgt, biases=bias))
#    it += 1
## Transfert Weights for Initial State
#NN.fit(X, y)
#NN.set_parameters(wgtandbias)
#check_wgt = NN.get_parameters()[0][0]
#check_wgt1 = wgtandbias[0][0]


## Predict With Saved State
np.random.seed(1)
#check_prob = NN.predict_proba(X)
## Reset
#reset = NN._backend._initialize_impl(X, y)
## Dimension of Multivariate Marginal Effects
n_inter = 2
## Create Containers for Multivariate Marginal Effects
cont_NN = create_out(n_inter)
cont_rf = create_out(n_inter)
## Compute Multivariate Effects Y = 1
c_mm = 0

## Bootstrap Marginal Effects w/ Standard Errors
for i in range(n_boot):
    boot(i, NN, X, y, marg_effects_NN, gof_measures_NN, "ETH Neural Network", cont_NN, n_inter, c_mm) #, weight=wgt)
    boot(i, rf, X, y, marg_effects_rf, gof_measures_rf, "ETH Random Forest", cont_rf, n_inter, c_mm)

out_NN = margeffects(marg_effects_NN, gof_measures_NN, n_boot, predictors)
out_rf = margeffects(marg_effects_rf, gof_measures_rf, n_boot, predictors)

## Output Variable Importance Random Forest
importances = rf.feature_importances_
std = np.std([tree.feature_importances_ for tree in rf.estimators_],
             axis=0)
indices = np.argsort(importances)[::-1]
import_features = np.concatenate((predictors[indices], importances[indices].reshape(len(features),1)), axis=1)

np.savetxt("./output/ETH_importance_rf.csv", import_features, fmt="%s", delimiter=',')
np.savetxt("./output/ETH_marginal_effects_NN.csv", out_NN, fmt="%s", delimiter=',')
np.savetxt("./output/ETH_GofF_NN.csv", gof_measures_NN, fmt="%s", delimiter=',')
np.savetxt("./output/ETH_marginal_effects_rf.csv", out_rf, fmt="%s", delimiter=',')
np.savetxt("./output/ETH_GofF_rf.csv", gof_measures_rf, fmt="%s", delimiter=',')

if c_mm == 1:
    out_MM_NN = margeffects_multi(cont_NN['out_by_x_m'], cont_NN['out_by_x_n'])
    out_MM_rf = margeffects_multi(cont_rf['out_by_x_m'], cont_rf['out_by_x_n'])
    np.savetxt("./output/ETH_marginal_mv_NN.csv", out_MM_NN, fmt="%s", delimiter=',')
    np.savetxt("./output/ETH_marginal_mv_rf.csv", out_MM_rf, fmt="%s", delimiter=',')



