###############################################################################
# svd.py
# Authors: Adam Fisch and Jonathan Shor
# COS 424 Assignment 3
#
# Description: Preliminary analysis of the SVD approach, prior to
# consolidating the code.
###############################################################################
from optparse import OptionParser
from GraphView import GraphView
import numpy as np
import Utils
import PerfMetrics
import time, sys, math
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt
from scipy.sparse.linalg import svds

###############################################################################
# A helper method to extract relevant data from the singular value
# decomposition approach. Returns the 30 largest singular values, their
# explained variances, and a 2 dimensional projection of the data
# (for easy investigative plotting).
###############################################################################
def svdData(Xtrain):
    (u, s, vt) = svd_svals = svds(Xtrain, k = 30)
    svd_variance = TruncatedSVD(n_components = 30)
    svd_variance.fit(Xtrain)
    variance_ratios = svd_variance.explained_variance_ratio_
    svd = TruncatedSVD(n_components = 2)
    svd.fit(Xtrain)
    X2D = svd.fit_transform(Xtrain)
    return (s, variance_ratios, X2D)

###############################################################################
# Helper method to run a single SVD decomposition, selecting only n
# eigenvectors. Returns the reconstructed, test vector.
###############################################################################
def doSVD(Xtrain, test, n):
    (u, s, vt) = svd_svals = svds(Xtrain, k = 30)
    test_points = test.shape[0]
    testPredict = [0]*test_points
    for i in range(test_points):
	    sender = test[i, 0]
	    receiver = test[i, 1]
	    predicted = 0
	    for dim in range(0, n):
		    u_component = u[sender, dim]
		    v_component = vt[dim, receiver]
		    predicted += (u_component * s[dim] * v_component)
	    testPredict[i] = predicted
    return (np.array(testPredict))

###############################################################################
# Run the steps for complete analysis of the SVD approach to a given sparse
# matrix Xsparse. Test is the matrix of probe data to validate on.
###############################################################################
def svdAnalysis(Xsparse, Test):
    # Visualize the distribution of the transaction data. Semilog scale in y
    # as the values are heavy tailed.
    plt.figure()
    (n, bins, patches) = plt.hist(Xsparse.data, 30, normed = 1)
    plt.yscale('log')

    # Compute the details on a SVD decomposition - the explained variance
    # ratios, the singular values themselves, and a 2D projection for plotting.
    (singular_vals, variance_ratios, X2d) = svdData(Xsparse)
    plt.figure()
    plt.plot(np.fliplr([singular_vals])[0], 'bo-')
    plt.title('Singular Values: Decreasing Magnitude')
    plt.figure()
    plt.plot(variance_ratios, 'bo-')
    plt.title('Singular Values: Decreasing Explained Variance')
    plt.figure()
    plt.scatter(X2d[:,0], X2d[:,1])
    plt.title('2D Projection of Transaction Data')

    fprs = []
    tprs = []
    aucs = []
    precisions = []
    recalls = []
    list_of_metrics = [fprs, tprs, aucs, precisions, recalls]

    # TODO: Cross validation for # of s vals.
    # Compute the actual low rank approximation through the selected SVD
    print ("doing regular method")
    Xpredict = doSVD(Xsparse, Test, 12)
    # Compute errors
    (fpr, tpr, auc, precision, recall) = \
        PerfMetrics.performance_metrics(Test, Xpredict)
    append_lists(list_of_metrics, [fpr, tpr, auc, precision, recall])

    print ("doing binary method")
    binaryXsparse = binarize(Xsparse);
    binaryXpredict = doSVD(binaryXsparse, Test, 12)
    # Compute errors
    (fpr, tpr, auc, precision, recall) = \
        PerfMetrics.performance_metrics(Test, binaryXpredict)
    append_lists(list_of_metrics, [fpr, tpr, auc, precision, recall])

    print ("doing bins method")
    binnedXsparse = bin_a_rize(Xsparse,12)
    binnedXpredict = doSVD(binnedXsparse, Test, 12)
    (fpr, tpr, auc, precision, recall) = \
        PerfMetrics.performance_metrics(Test, binnedXpredict)
    append_lists(list_of_metrics, [fpr, tpr, auc, precision, recall])

    PerfMetrics.plotROC(list_of_metrics[0], list_of_metrics[1])
    PerfMetrics.plotPR(list_of_metrics[3], list_of_metrics[4])

    print aucs
    return

###############################################################################
# Update list of metrics.
###############################################################################
def append_lists(list_lists, list_vals):
    for i in range(len(list_lists)):
        list_lists[i].append(list_vals[i])

    return list_lists

###############################################################################
# Count transactions as either true or false.
###############################################################################
def binarize(Xsparse):
    X_binary = sp.coo_matrix(Xsparse)
    X_binary.data[:] = 1
    return (X_binary)

###############################################################################
# Bin transactions counts. Put it in a bin according to what power of e it is.
# This should make it a bit denser near the bottom. Add 1, so that you start
# bin counts at 1. 0 is the auto bin for the sparse entries.
###############################################################################
def bin_a_rize(Xsparse, n):
    X_bins = sp.coo_matrix(Xsparse)
    values = X_bins.data
    for i in range(values.shape[0]):
        bin = math.ceil(math.log(values[i])) + 1
        if (bin > n):
            bin = n
        values[i] = bin
    return (X_bins)

###############################################################################
# Main. Read in data and run analysis.
###############################################################################
def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path
    start_time = time.time()

    Test = Utils.read_test_trps_txt(path)
    Xtrain = Utils.read_train_trps_txt(path)
    Xtrain_sparse = GraphView.get_coo(Xtrain)
    print Xtrain_sparse
    svdAnalysis(Xtrain_sparse, Test)
    print time.time() - start_time
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])
