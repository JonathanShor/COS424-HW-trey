'''
Created on Apr 13, 2015

@author: jonathanshor
'''
from optparse import OptionParser
# from Explore import Explore as ex
import time,sys
import numpy as np
import networkx as nx
import scipy.sparse as sp

# TRAIN_FNAME = "txTripletsCountsNo_e-05.txt"
TRAIN_FNAME = "txTripletsCountsWGiantOnly.txt"
TEST_FNAME = "testTriplets.txt"

def read_train_trps_txt(path, toNX=False, skip = 0, fn = ''):
# Accepts path to TRAIN_FNAME
# If toNX, returns a nx.DiGraph, otherwise returns a ndarray
# Can be given a number of rows to skip, ndarray case only
    if len(fn)>0:
        fname = fn
    else:
        fname = TRAIN_FNAME
    if toNX:
        return nx.read_weighted_edgelist(path + fname, create_using=nx.DiGraph(), nodetype=int)
    return np.loadtxt(path + fname, skiprows = skip)

def read_test_trps_txt(path, toNX=False, skip = 0):
# Accepts path to TEST_FNAME
# If toNX, returns a nx.DiGraph, otherwise returns a ndarray
# Can be given a number of rows to skip, ndarray case only
    if toNX:
        return nx.read_weighted_edgelist(path + TEST_FNAME, create_using=nx.DiGraph(), nodetype=int)
    return np.loadtxt(path + TEST_FNAME, skiprows = skip)

def get_coo(trp_raw):
# Given nX3 raw trp ndarray trp_raw, return a sparse.coo_matrix of minimal square shape
    dim = max(trp_raw[:,0].max(),trp_raw[:,1].max()) + 1
    return sp.coo_matrix((trp_raw[:,2],(trp_raw[:,0],trp_raw[:,1])),shape=(dim,dim))

def dedup(raw):
# Takes raw Nx3 triplet format matrix, checks for duplicate sender-receiver entries
# Returns Mx3 ndarray, with the dup entries removed
# NOT TESTED
    keys = list("+".join([str(raw[i][0]),str(raw[i][1])]) for i in range(len(raw)))
    return raw[np.unique(keys,return_index=True)[1]]

def make_held_out(train, alph=0.1):
# Given a training set (nX3 ndarray of trpl) and percent alph
# Returns disjoint subsets of train, splitting train into size alph% and (1-alph)% subsets
# Size is taken to be on number of transactions
    print "train size: %s, alph=%s" % (train[:,2].sum(),alph)
    trainset = train.copy()
    valset = np.empty((0,3),dtype='int32')
    for i in range(len(trainset)):
        delt = trainset[i,2] * alph
        delt = int(delt) + (delt % 1 > np.random.rand())
        if delt:
            trainset[i,2] -= delt
            valset = np.append(valset, [[trainset[i,0],trainset[i,1],delt]],0)
    trainset = trainset[trainset[:,2] != 0]
    print "trainset: %s  valset: %s total: %s" % (trainset[:,2].sum(), valset[:,2].sum(), trainset[:,2].sum() + valset[:,2].sum())
    return (valset, trainset)

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read data from PATH', metavar='PATH')
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path

    start_time = time.time()

#     test = read_test_trps_txt(path)
#     print "Test.shape = %s" % str(test.shape)
    if options.path:
        train = read_train_trps_txt(path)
        print "Train.shape = %s" % str(train.shape)

#     ex.collect_alt_views(dedup(test), path + 'DEDUPtestTriplets.txt', \
#                               comments='Duplicates entries removed.')
#     ex.collect_alt_views(dedup(train), path + 'DEDUPtxTripletsCounts.txt', \
#                               comments='Duplicates entries removed.')


    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])
