'''
Created on Apr 20, 2015

@author: jonathanshor
'''
import time, sys
from optparse import OptionParser
import Utils
import scipy.sparse.csgraph

def get_coo(trp_raw):
# Given nX3 raw trp ndarray trp_raw, return a sparse.coo_matrix of minimal square shape
    dim = max(trp_raw[:,0].max(),trp_raw[:,1].max()) + 1
    return scipy.sparse.coo_matrix((trp_raw[:,2],(trp_raw[:,0],trp_raw[:,1])),shape=(dim,dim))

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path

    start_time = time.time()

    train = get_coo(Utils.read_train_trps_txt(path))
    print "Connected comps follow"
    print scipy.sparse.csgraph.connected_components(train, directed=True, connection='strong', return_labels=True)
    
    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])