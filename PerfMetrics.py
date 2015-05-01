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
#       matrix: The predicted matrix of transactions
#
# Output: FPR, TPR, AUC, F1, P, R
###############################################################################
def performance_metrics(test_entries, matrix):
    # Gather and process data
    matrix_prediction = constrain(matrix)
    test_labels = test_entries[:,2]
    predictions = []
    for i in range(test_labels.shape[0]):
        sender = test_entries[i,0]
        receiver = test_entries[i,1]
        predictions[i] = matrix_prediction[sender, receiver]
    predictions = np.array(predictions)

    # roc_curve returns false-positive and true-positive rates
    (fpr, tpr, _) = roc_curve(test_labels, predictions)

    # AUC (area under ROC curve)
    auc = auc(test_labels, predictions)

    # f1 (combination of precision and recall)
    f1 = f1_score(test_labels, predictions)

    # precision/recall estimates
    (precision, recall, _) = precision_recall_curve(test_labels, predictions)

    return (fpr, tpr, auc, f1, precision, recall)

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
def constrain(matrix):
    rows = matrix.shape[0]
    cols = matrix.shape[1]
    for i in range(rows):
        for j in range(cols):
            if (matrix[i, j] > 1):
                matrix[i, j] = 1
            elif (matrix[i, j] < 0):
                matrix[i, j] = 0

    return (matrix)

###############################################################################
# Helper method to plot the results of the ROC curve.
###############################################################################
def plotROC(fpr, tpr):
    ax = plt.figure()
    colors = ['b','g','r','c','m','y','k']
    for c in range(len(fpr)):
        (this_fpr, this_tpr) = (fpr[c], tpr[c])
        plt.plot(this_fpr, this_tpr, colors[c] + '-')

    x = np.linspace(0,1,num = len(fpr))
    plt.plot(x, x, '--', color = (.6, .6, .6), label = 'Luck')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Reciever Operating Characteristic')
    plt.legend(loc = 'best', prop={'size':10})

    return



