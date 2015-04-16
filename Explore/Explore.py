'''
Created on Apr 13, 2015

@author: jonathanshor
'''
from optparse import OptionParser
import time, sys
import numpy as np
import Utils
import matplotlib.pyplot as plt

def stats_explore(vals):
# Requires a nX1 ndarray
    print "Mean: %0.2f, standard dev: %0.2f" % (np.mean(vals), vals.std())
    print "Median: %0.2f" % np.median(vals)
    x = np.linspace(vals.min(),vals.max(),vals.max())
    print "Range: [%s, %s]" % (x[0], x[-1])
    
    y = x.copy()
    y = [len(vals[np.where(vals == i)]) for i in y]
    print "Non-zero ys: %s" % np.count_nonzero(y)
    plt.plot(x,y,'o')
    plt.show()

def trans_freq_1way_view(trp):
# Requires nx3 ndarray from a triplets file
# Returns a mX2 ndarray:
# column 1 contains the count of how many address pairs had exactly the number in column 0 one-way transactions
    trp = trp[:,2]
    x = np.linspace(trp.min(),trp.max(),trp.max())
    y = [len(trp[np.where(trp == i)]) for i in x]
    view = np.empty((len(x),2))
    view[:,0] = x
    view[:,1] = y
    view = view[view[:,1].nonzero()]
    return view

#def trans_freq_2way_view(trp):
# Requires nx3 ndarray from a triplets file
# Returns a mX2 ndarray:
# column 1 contains the count of how many address pairs had exactly the number in column 0 transactions
# This combines transactions in both directions


def collect_alt_views(array, name, comments = ""):
# Requires an ndarray and a filename (with path) and optional comments for the file
    np.savetxt(name, array, header = comments)

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path

    start_time = time.time()

    train = Utils.read_train_trps_txt(path)
    #print "Transaction stats"
    #stats_explore(train[:,2])
    collect_alt_views(trans_freq_1way_view(train), path + 'transFreqView.txt', \
                       'Number of transactions; Count of one way address pairs')
    

    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])