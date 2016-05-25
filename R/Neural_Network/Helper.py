# -*- coding: utf-8 -*-
# Helper Functions
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# May 25, 2016

from sknn.mlp import Classifier, Layer
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.grid_search import GridSearchCV
from sklearn import cross_validation
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
import numpy as np

########################################################################      
## ROC - AUC                                                           #
########################################################################
def scorer(estimator, xt, yt, xte, yte):
    # Train the Classifier
    # estimator.fit(xt, yt)
    # Predict Probabilities
    np.random.seed(1)
    Prob = estimator.predict_proba(xte)    
    # Compute ROC and AUC
    ROC = roc_curve(yte, Prob[:, 1])
    fpr = ROC[0]
    tpr = ROC[1]
    thr = ROC[2]
    # Extract Threshold Maximizing Sensitivity + Specifity
    crit = np.argmax((1-fpr) + tpr)
    max_thr = thr[crit]
    AUC = auc(fpr, tpr)
    return AUC, fpr, tpr, max_thr

########################################################################      
# Plot of a ROC curv                                                   #
########################################################################
def pics(auc, name, rep):
    with PdfPages('%s' %name+str(rep) + '.pdf') as pdf:
        plt.figure()
        plt.plot(auc[1], auc[2], label='ROC curve (area = %0.3f)' % auc[0])
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic %s' % name )  
        plt.legend(loc="lower right")
        #plt.show()
        pdf.savefig()
        plt.close

########################################################################
# Marginal Effets                                                      #
########################################################################
def marg(estimator, boot, smpl, out):
    n_effect = len(features)
    for j in range(n_effect):
        xtest0 = smpl.copy()
        xtest0[:,j] = xtest0[:,j] * 0.0
        xtest1 = smpl.copy()
        xtest1[:,j] = xtest1[:,j] * 0.0 + 1.0
        part1 = np.mean(estimator.predict_proba(xtest1)[:,1])
        part2 = np.mean(estimator.predict_proba(xtest0)[:,1])
        out[j,boot] =  (part1 - part2)

########################################################################
# Meta Parameters Calibration function                                 #
########################################################################
def meta(X_, y_, act):
    ## Neural Network Classifier -- 3 Hidden Layer
    NN = Classifier(layers = [Layer(act, units=20), 
                              Layer(act, units=20),
                              Layer("Softmax")],
                              regularize="L2",
                              n_iter = 1000,
                              verbose=True,
                              batch_size=32,
                              learning_rule="adagrad",
                              random_state=0)
    ## Meta Parameters Grid Search with Cross Validation
    param_grid = {"learning_rate": [0.001, 0.01],
                  "weight_decay": [0.0001, 0.001],
                  "hidden0__units": [50, 100],
                  "hidden1__units": [50, 100]}

    NN = GridSearchCV(NN, param_grid, refit=True, verbose=True, scoring='roc_auc', n_jobs=1, cv=5)
    ## Fit the Classifier
    np.random.seed(1)
    NN.fit(np.asarray(X_), np.asarray(y_, dtype=np.int8))
    ## Best Fit Estimator
    Best_NN = NN.best_estimator_
    
    ## Random Forest
    rf = RandomForestClassifier(random_state=0, verbose=1, n_estimators=1000)
    ## Meta Parameters Grid Search with Cross Validation
    param_grid = {'max_features': ["auto", "log2", np.int(np.shape(X_)[1]/2)],
                  'n_estimators': [100,500,1000]}    
    rf = GridSearchCV(rf, param_grid, refit=True, verbose=True, scoring='roc_auc', n_jobs=1, cv=5)
    ## Fit the Classifier
    np.random.seed(1)
    rf.fit(np.asarray(X_), np.asarray(y_, dtype=np.int8))
    ## Best Fit Estimator
    Best_rf = rf.best_estimator_
    #Best_NN,
    return [Best_NN, Best_rf]

########################################################################
# Auto-Encoder                                                         #
########################################################################
def auto(X_, act, units_):
    ## Neural Network Classifier -- 3 Hidden Layer
    myae = Classifier(layers = [Layer(act, units=units_[0]), 
                                Layer("Softmax")],
                                n_iter = 100,
                                verbose=True,
                                regularize="L2",
                                batch_size=32,
                                learning_rule="adagrad")
    ## Fit the Classifier
    np.random.seed(1)
    myae.fit(np.asarray(X_), np.asarray(X_))
    
    return myae


########################################################################
# Bootstrap function                                                   #
########################################################################
def boot(rep, estimator, X_, y_, out, gof, clas, wgt_):

    ## Split Train Test for Marginal Effects Bootstrap
    X_train, X_test, y_train, y_test = cross_validation.train_test_split(np.asarray(X_), np.asarray(y_), test_size=0.33, random_state=rep)
    ## AUC and Sensitivity on Left Out Test 
    np.random.seed(1)
    ## Transfer Weights
    #if estimator == "NN":
    #    estimator.set_parameters(wgt_)
    estimator.fit(X_train, y_train)
    AUC = scorer(estimator, X_train, y_train, X_test, y_test)
    ## Test Set Predictions -- Sensitivity and Specificity
    np.random.seed(1)
    y_prob_pred = estimator.predict_proba(X_test)[:,1].reshape(len(y_test),1)
    # Categorize Predictions
    y_pred = np.zeros((len(y_test),1))
    y_pred[:,0]  = y_prob_pred[:,0] >= AUC[3]
    y_tt = y_test.reshape(len(y_test),1)
    NP = np.sum(y_tt)
    NG = np.sum(1 - y_tt)
    TP = np.sum(np.multiply(y_tt, y_pred))
    TN = np.sum(np.multiply((1-y_tt), (1-y_pred)))
    sens = np.float(TP) / np.float(NP)
    spec = np.float(TN) / np.float(NG)

    gof[0,rep] = AUC[0]
    gof[1,rep] = sens
    gof[2,rep] = spec
    gof[3,rep] = AUC[3]

    ## Compute Marinal Effect
    marg(estimator, rep, X_test, out)
    ## Export AUC chart every 10 bootstrap
    if rep % 10 == 0:
        pics(AUC, clas, rep)

########################################################################
# Marginal Effects                                                     #
########################################################################
def margeffects(marg_effects, gof_measures, boot, predictors):
    marg_effects[:, boot] = np.mean(marg_effects[:,range(boot)], axis=1) 
    marg_effects[:, (boot+1)] = np.std(marg_effects[:,range(boot)], axis=1) 
    marg_effects[:, (boot+2)] = marg_effects[:, boot] / marg_effects[:, (boot+1)]
    gof_measures[:, boot] = np.mean(gof_measures[:,range(boot)], axis=1)
    gof_measures[:, (boot+1)] = np.std(gof_measures[:,range(boot)], axis=1)
    ## Export Results
    toout = np.concatenate([predictors, marg_effects], axis=1)
    ind = np.abs(np.double(toout[:, (boot+3)])).argsort()[::-1]
    out = toout[ind]
    return out
