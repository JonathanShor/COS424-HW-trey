'''
Created on Apr 13, 2015

@author: jonathanshor
'''
from optparse import OptionParser
import time,sys
import numpy as np
import networkx as nx

TRAIN_FNAME = "txTripletsCounts.txt"
TEST_FNAME = "testTriplets.txt"

def read_train_trps_txt(path, toNX=False, skip = 0):
# Accepts path to TRAIN_FNAME
# If toNX, returns a nx.DiGraph, otherwise returns a ndarray
# Can be given a number of rows to skip, ndarray case only
    if toNX:
        return nx.read_weighted_edgelist(path + TRAIN_FNAME, create_using=nx.DiGraph(), nodetype=float)
    return np.loadtxt(path + TRAIN_FNAME, skiprows = skip)

def read_test_trps_txt(path, toNX=False, skip = 0):
# Accepts path to TEST_FNAME
# If toNX, returns a nx.DiGraph, otherwise returns a ndarray
# Can be given a number of rows to skip, ndarray case only
    if toNX:
        return nx.read_weighted_edgelist(path + TRAIN_FNAME, create_using=nx.DiGraph(), nodetype=float)
    return np.loadtxt(path + TEST_FNAME, skiprows = skip)

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

    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])