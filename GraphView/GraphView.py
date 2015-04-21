'''
Created on Apr 20, 2015

@author: jonathanshor
'''
import time, sys
from optparse import OptionParser
import Utils
import Explore.Explore as ex
import scipy.sparse.csgraph as csg
import scipy.sparse as sp
#import networkx as nx

def get_coo(trp_raw):
# Given nX3 raw trp ndarray trp_raw, return a sparse.coo_matrix of minimal square shape
    dim = max(trp_raw[:,0].max(),trp_raw[:,1].max()) + 1
    return sp.coo_matrix((trp_raw[:,2],(trp_raw[:,0],trp_raw[:,1])),shape=(dim,dim))

# def get_coSCC_matrix(comp_labs, full_label=False):
# Requires comp_labs = len N list of SCC labels (expected from scipy.sparse.csgraph.connected_components)
# Returns an NxN scipy.sparse.lil_matrix such that position i,j is 1 iff vertices /
# i and j are in the same SCC and 0 otherwise.
# If full_label=True, position i,j contains the root (labeling) vertex of the SCC + 1, instead of 1.
# Note the +1 to the label value; labeling starts at 0. An entry of 0 still means not in the same SCC. 
    #co_ssc = sp.lil_matrix((len(comp_labs),len(comp_labs)))
#     co_ssc = sp.dok_matrix((len(comp_labs),len(comp_labs)))
#     for i in range(len(comp_labs)):
#         pass

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
    parser.add_option("-s", dest="strongcc", action="store_true", default=False)
    parser.add_option("-w", dest="weakcc", action="store_true", default=False)
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path

    start_time = time.time()

    train = get_coo(Utils.read_train_trps_txt(path))
    
    if options.weakcc:
        (_numSCCs, scc_labels) = csg.connected_components(train, directed=True, \
                                                      connection='weak', return_labels=True)
        ex.collect_alt_views(ex.SCCs_freq_view(scc_labels), path + "WCCsCountView.txt", \
                         comments="WCC label (arbitrary); Count of vertex in WCC")
    if options.strongcc:
        (_numSCCs, scc_labels) = csg.connected_components(train, directed=True, \
                                                      connection='strong', return_labels=True)
        ex.collect_alt_views(ex.SCCs_freq_view(scc_labels), path + "SCCsCountView.txt", \
                         comments="SCC label (arbitrary); Count of vertex in SCC")
    
    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])