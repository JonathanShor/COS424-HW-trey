###############################################################################
# PerfMetrics.py
# Authors: Adam Fisch and Jonathan Shor
# COS 424 Assignment 3
#
# Description: Calculate performance metrics for our predictions.
###############################################################################
from sklearn.metrics import *
import scipy.sparse as sp
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Calculate performance metrics.
# Input:
#       test_entries: The sparse matrix of testing probes
#       test_predictions: The predicted matrix of testing probes
#
# Output: FPR, TPR, AUC, F1, P, R
###############################################################################
def performance_metrics(test_entries, test_predictions):
    # Gather and process data
    predictions = constrain(test_predictions)
    test_labels = test_entries[:,2]

    # roc_curve returns false-positive and true-positive rates
    (fpr, tpr, _) = roc_curve(test_labels, predictions)

    # AUC (area under ROC curve)
    aucScore = auc(fpr, tpr)

    # precision/recall estimates
    (precision, recall, _) = precision_recall_curve(test_labels, predictions)

    return (fpr, tpr, aucScore, precision, recall)

###############################################################################
# Helper method to change values to weighted probabilities of sender i giving
# payment to receiver j.
###############################################################################
def makeWeightedProbabilities(matrix):
    rows = matrix.shape[0]
    cols = matrix.shape[1]
    col_sums = matrix.sum(1)
    for i in range(rows):
        sum = col[i,0]
        if (sum == 0):
            continue
        for j in range(cols):
            matrix[i, j] = matrix[i, j]/sum

    return (matrix)

###############################################################################
# Helper method that just take the real value of the matrix. Any transaction
# prediction above 1 is taken as truth. Any negative values are treated as
# absolute no (0). And in between (0,1) are probabilities.
###############################################################################
def constrain(array):
    for i in range(len(array)):
        if (array[i] > 1):
            array[i] = 1
        elif (array[i] < 0):
            array[i] = 0

    return (array)

###############################################################################
# Helper method to plot the results of the ROC curve.
###############################################################################
def plotROC(fpr, tpr):
    #print fpr.shape
    print len(fpr)
    ax = plt.figure()
    colors = ['b','g','r','c','m','y','k']
    labels = ['Standard', 'Binary', 'Binned']
    for c in range(len(fpr)):
        (this_fpr, this_tpr) = (fpr[c], tpr[c])
        plt.plot(this_fpr, this_tpr, colors[c] + '-', label = labels[c])

    x = np.linspace(0,1,num = len(fpr))
    plt.plot(x, x, '--', color = (.6, .6, .6), label = 'Luck')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Reciever Operating Characteristic')
    plt.legend(loc = 'best', prop={'size':10})

    return

###############################################################################
# Helper method to plot the results of the PR curve.
###############################################################################
def plotPR(precision, recall):
    #print fpr.shape
    print len(precision)
    ax = plt.figure()
    colors = ['b','g','r','c','m','y','k']
    labels = ['Standard', 'Binary', 'Binned']
    for c in range(len(precision)):
        (this_pr, this_r) = (precision[c], recall[c])
        plt.plot(this_pr, this_r, colors[c] + '-', label = labels[c])

    x = np.linspace(0,1,num = len(precision))
    plt.plot(x, x, '--', color = (.6, .6, .6), label = 'Luck')
    plt.xlabel('Precision')
    plt.ylabel('Recall')
    plt.title('Precision - Recall Curve')
    plt.legend(loc = 'best', prop={'size':10})

    return



