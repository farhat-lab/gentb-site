# -*- coding: utf-8 -*-
# Helper Functions
# Author: Jimmy Royer
# jimmy.royer@analysisgroup.com
# May 15, 2016

from sknn.mlp import Classifier, Layer
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.grid_search import GridSearchCV
from sklearn import cross_validation
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import recall_score

########################################################################      
## ROC - AUC                                                           #
########################################################################
def scorer(estimator, xt, yt, xte, yte):
    # Train the Classifier
    estimator.fit(xt, yt)
    # Predict Probabilities
    Prob = estimator.predict_proba(xte)    
    # Compute ROC and AUC
    ROC = roc_curve(yte, Prob[:, 1])
    fpr = ROC[0]
    tpr = ROC[1]
    AUC = auc(fpr, tpr)
    return AUC, fpr, tpr
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
def marg(estimator, boot, smpl):
    n_effect = len(features)
    for j in range(n_effect):
        xtest0 = smpl.copy()
        xtest0[:,j] = xtest0[:,j] * 0.0
        xtest1 = smpl.copy()
        xtest1[:,j] = xtest1[:,j] * 0.0 + 1.0
        part1 = np.mean(estimator.predict_proba(xtest1)[:,1])
        part2 = np.mean(estimator.predict_proba(xtest0)[:,1])
        marg_effects[j,boot] =  (part1 - part2)
########################################################################
# Bootstrap function                                                   #
########################################################################
def boot(rep):
    
    ## Split Train Test for Marginal Effects Bootstrap
    X_train, X_test, y_train, y_test = cross_validation.train_test_split(np.asarray(X), np.asarray(y), test_size=0.33, random_state=rep)

    ## Neural Network Classifier -- 2 Hidden Layer
    NN = Classifier(layers = [Layer("Rectifier", units=20), 
                              Layer("Rectifier", units=20), 
                              Layer("Softmax")],
                              regularize="L2",
                              n_iter = 1000,
                              verbose=True,
                              batch_size=32,
                              learning_rule="adagrad",
                              random_state=12346)

    ## Meta Parameters Grid Search with Cross Validation
    param_grid = {"learning_rate": [0.001, 0.01],
                  "weight_decay": [0.0001, 0.001],
                  "hidden0__units": [25, 50, 100],
                  "hidden1__units": [25, 50, 100]}
    NN = GridSearchCV(NN, param_grid, refit=True, verbose=True, scoring='roc_auc', n_jobs=1, cv=3)

    ## Fit the Classifier
    np.random.seed(1)
    NN.fit(X_train, np.asarray(y_train, dtype=np.int8))

    ## Best Fit Estimator
    Best = NN.best_estimator_

    ## Test Set Predictions
    y_pred = Best.predict(X_test)

    ## AUC and Sensitivity on Left Out Test 
    np.random.seed(1)
    AUC = scorer(Best, X_train, y_train, X_test, y_test)
    sensitivity = recall_score(y_test, y_pred)
    gof_measures[0,rep] = AUC[0]
    gof_measures[1,rep] = sensitivity

    ## Compute Marinal Effect
    np.random.seed(1)
    marg(Best, rep, X_test)
    
    ## Export AUC chart every 10 bootstrap
    if rep % 10 == 0:
        pics(AUC, "Neural_Network",rep)
