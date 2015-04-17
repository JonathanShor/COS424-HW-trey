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
    
def sender_rank_view(trp, receivers=False ):
# Requires nx3 ndarray from a triplets file -- MAY BE MODIFIED
# Optional receivers bool to return receivers rank view
# Returns a mX2 ndarray: Address and (sorted by) number of sends(receives) 
    if receivers:
        trp = trp.sort(axis=1)
        addrs = trp[:,1]
    else:
        # Assumes sorted by senders already 
        addrs = trp[:,0]
    senders = np.unique(addrs)
#    sends = [trp[:,2][np.where(addrs == x)].sum() for x in senders]
    sends = [trp[:,2][np.searchsorted(trp[:,0], x, 'left'): \
                      np.searchsorted(trp[:,0], x, 'right')].sum() for x in senders]
    print "Senders: %s, sends: %s" % (len(senders), len(sends))
    ret = np.empty((len(senders),2))
    ret[:,0] = senders
    ret[:,1] = sends
#    return ret.sort(kind='mergesort')
    return ret

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
    collect_alt_views(sender_rank_view(train.copy()), path + "SendersRanked.txt", \
                      comments="Sender Address; Count of Sends")
#     collect_alt_views(sender_rank_view(train.copy(), receivers=True), path + "ReceiversRanked.txt", \
#                       comments="Receiver Address; Count of Receipts")

    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])