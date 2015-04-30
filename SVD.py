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
import time, sys
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
# eigenvectors. Returns the reconstructed, low-rank matrix.
###############################################################################
def doSVD(Xtrain, n):
    svd = TruncatedSVD(n_components = n)
    svd.fit(Xtrain)
    Xpredict = svd.fit_transform(Xtrain)
    return (Xpredict)

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
    print singular_vals.shape
    plt.plot(np.fliplr([singular_vals])[0], 'bo-')
    plt.title('Singular Values: Decreasing Magnitude')
    plt.figure()
    plt.plot(variance_ratios, 'bo-')
    plt.title('Singular Values: Decreasing Explained Variance')
    plt.figure()
    plt.scatter(X2d[:,0], X2d[:,1])
    plt.title('2D Projection of Transaction Data')

    # TODO: Cross validation for # of s vals.
    # Compute the actual low rank approximation through the selected SVD
    Xpredict = doSVD(Xsparse, 12)

    # Compute errors
    # TODO: Threshold for a valid transaction?
    errors = None
    return errors

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

    Xtrain = Utils.read_train_trps_txt(path)
    Xtrain_sparse = GraphView.get_coo(Xtrain)

    svdAnalysis(Xtrain_sparse, None)
    plt.show()



if __name__ == '__main__':
    main(sys.argv[1:])