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
    trans = trp[:,2]
    x = np.linspace(trans.min(),trans.max(),trans.max())
    y = [len(trans[np.where(trans == i)]) for i in x]
    view = np.empty((len(x),2))
    view[:,0] = x
    view[:,1] = y
    view = view[view[:,1].nonzero()]
    return view

def SCCs_freq_view(SCCs):
# Given SCCs labels array of len N, return Nx2 ndarray with counts of each SCC label
# Sorted
    num_SCC = len(np.unique(SCCs))
    view = np.empty((num_SCC,2))
    view[:,0] = np.arange(num_SCC)
    view[:,1] = [len(SCCs[np.where(SCCs == i)]) for i in view[:,0]]
    return view[np.argsort(view[:,1])][::-1]

# def trans_freq_2way_view(trp):
# Requires nx3 ndarray from a triplets file
# Returns a mX2 ndarray:
# column 1 contains the count of how many address pairs had exactly the number in column 0 transactions
# This combines transactions in both directions
#     sends = trp.copy()
#     recs = trp
#     recs[:,0] = 
#     for i in range(len(trp)):
#          
#     trans = trp[:,2]
#     x = np.linspace(trans.min(),trans.max(),trans.max())
    
def sender_rank_view(trips, receivers=False ):
# Requires nx3 ndarray from a triplets file -- MAY BE MODIFIED
# Optional receivers bool to return receivers rank view instead of default senders
# Returns a mX2 ndarray: Address and (rev sorted by) number of sends(receives)
# m = # unique senders (receivers)
    if receivers:
        # Sort by senders, and swap senders into column 0, rec into col 1
        temp = trips[np.argsort(trips[:,1])]
        trp = temp.copy()
        trp[:,0] = temp[:,1]
        trp[:,1] = temp[:,0]
#        addrs = trp[:,1]
    else:
        # Assumes sorted by senders already
        trp = trips 
    addrs = trp[:,0]
    senders = np.unique(addrs)
#    sends = [trp[:,2][np.where(addrs == x)].sum() for x in senders]
    sends = [trp[:,2][np.searchsorted(trp[:,0], x, 'left'): \
                      np.searchsorted(trp[:,0], x, 'right')].sum() for x in senders]
#    print "Senders: %s, sends: %s" % (len(senders), len(sends))
#    print "Sum matching, orig -- new: %s -- %s" % (trips[:,2].sum(), np.array(sends).sum())
    ret = np.empty((len(senders),2))
    ret[:,0] = senders
    ret[:,1] = sends
    return ret[np.argsort(ret[:,1])][::-1]

def CCgen_view(ccgen):
# Given networkx CC generator, return # of CCsX2 ndarray
# Col 0: one vertex in CC; Col 1 (sorted): # of vertex in CC
    first = ccgen.next()
    cc_cnts = np.array([[first[0],len(first)]])
    for cc in ccgen:
        cc_cnts = np.vstack((cc_cnts, [[cc[0], len(cc)]]))
    return cc_cnts[np.argsort(cc_cnts[:,1])][::-1]

def collect_alt_views(array, name, comments = ""):
# Requires an ndarray and a filename (with path) and optional comments for the file
    np.savetxt(name, array,fmt='%d', header = comments)

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
#     collect_alt_views(trans_freq_1way_view(train.copy()), path + 'trans1WayFreqView.txt', \
#                        'Number of transactions; Count of one way address pairs')
    print "Collecting SendersRankedView.txt"
    collect_alt_views(sender_rank_view(train.copy()), path + "SendersRankedView.txt", \
                      comments="Sender Address; Count of Sends")
    print "Collecting ReceiversRankedView.txt"
    collect_alt_views(sender_rank_view(train.copy(), receivers=True), path + "ReceiversRankedView.txt", \
                       comments="Receiver Address; Count of Receipts")

    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])