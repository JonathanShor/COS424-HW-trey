'''
Created on Apr 20, 2015

@author: jonathanshor
'''
import time, sys
from optparse import OptionParser
import Utils
from Explore import Explore as ex
import scipy.sparse.csgraph as csg
import scipy.sparse as sp
import numpy as np
import networkx as nx

def get_coo(trp_raw):
# Given nX3 raw trp ndarray trp_raw, return a sparse.coo_matrix of minimal square shape
    dim = max(trp_raw[:,0].max(),trp_raw[:,1].max()) + 1
    return sp.coo_matrix((trp_raw[:,2],(trp_raw[:,0],trp_raw[:,1])),shape=(dim,dim))

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
    parser.add_option("-s", dest="strongcc", action="store_true", default=False)
    parser.add_option("-w", dest="weakcc", action="store_true", default=False)
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path

    start_time = time.time()

#     train = get_coo(Utils.read_train_trps_txt(path))
    trainx = Utils.read_train_trps_txt(path, toNX=True)
    print "%.2f -- nx Graph is of len: %s" % (time.time()-start_time, len(trainx))
    
    if options.weakcc:
        print "%.2f -- Collecting WCCs" % (time.time()-start_time)
#         (_numSCCs, scc_labels) = csg.connected_components(train, directed=True, \
#                                                       connection='weak', return_labels=True)
#         ex.collect_alt_views(ex.SCCs_freq_view(scc_labels), path + "WCCsCountView.txt", \
#                          comments="WCC label (arbitrary); Count of vertex in WCC")
        cc_gen = nx.weakly_connected_components(trainx)
        ex.collect_alt_views(ex.CCgen_view(cc_gen), path + "WCCsXCountView.txt", \
                             comments= "Vertex from WCC; Count of vertex in WCC")
    if options.strongcc:
        print "%.2f -- Collecting SCCs" % (time.time()-start_time)
        cc_gen = nx.strongly_connected_components(trainx)
        ex.collect_alt_views(ex.CCgen_view(cc_gen), path + "SCCsXCountView.txt", \
                             comments= "Vertex from SCC; Count of vertex in SCC")
#         (_numSCCs, scc_labels) = csg.connected_components(train, directed=True, \
#                                                       connection='strong', return_labels=True)
#         ex.collect_alt_views(ex.SCCs_freq_view(scc_labels), path + "SCCsCountView.txt", \
#                          comments="SCC label (arbitrary); Count of vertex in SCC")
    
    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])