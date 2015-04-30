'''
Created on Apr 13, 2015

@author: jonathanshor
'''
from optparse import OptionParser
# from Explore import Explore as ex
import time,sys
import numpy as np
import networkx as nx

# TRAIN_FNAME = "txTripletsCountsNo_e-05.txt"
TRAIN_FNAME = "txTripletsCountsWGiantOnly.txt"
TEST_FNAME = "testTriplets.txt"

def read_train_trps_txt(path, toNX=False, skip = 0):
# Accepts path to TRAIN_FNAME
# If toNX, returns a nx.DiGraph, otherwise returns a ndarray
# Can be given a number of rows to skip, ndarray case only
    if toNX:
        return nx.read_weighted_edgelist(path + TRAIN_FNAME, create_using=nx.DiGraph(), nodetype=int)
    return np.loadtxt(path + TRAIN_FNAME, dtype='int32', skiprows = skip)

def read_test_trps_txt(path, toNX=False, skip = 0):
# Accepts path to TEST_FNAME
# If toNX, returns a nx.DiGraph, otherwise returns a ndarray
# Can be given a number of rows to skip, ndarray case only
    if toNX:
        return nx.read_weighted_edgelist(path + TEST_FNAME, create_using=nx.DiGraph(), nodetype=int)
    return np.loadtxt(path + TEST_FNAME, dtype='int32', skiprows = skip)

def dedup(raw):
# Takes raw Nx3 triplet format matrix, checks for duplicate sender-receiver entries
# Returns Mx3 ndarray, with the dup entries removed
# NOT TESTED
    keys = list("+".join([str(raw[i][0]),str(raw[i][1])]) for i in range(len(raw)))
    return raw[np.unique(keys,return_index=True)[1]]
    
    
def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read data from PATH', metavar='PATH')
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path

    start_time = time.time()

    test = read_test_trps_txt(path)
    print "Test.shape = %s" % str(test.shape)
    train = read_train_trps_txt(path)
    print "Train.shape = %s" % str(train.shape)
    
#     ex.collect_alt_views(dedup(test), path + 'DEDUPtestTriplets.txt', \
#                               comments='Duplicates entries removed.')
#     ex.collect_alt_views(dedup(train), path + 'DEDUPtxTripletsCounts.txt', \
#                               comments='Duplicates entries removed.')


    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])