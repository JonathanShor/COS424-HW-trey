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
# Take a nX1 ndarray
    print "Mean: %0.2f, standard dev: %0.2f" % (np.mean(vals), vals.std())
    print "Median: %0.2f" % np.median(vals)
    x = np.linspace(vals.min(),vals.max(),vals.max())
    print "Range: [%s, %s]" % (x[0], x[-1])
    
    y = x.copy()
    y = [len(vals[np.where(vals == i)]) for i in y]
    print "Non-zero ys: %s" % np.count_nonzero(y)
    plt.plot(x,y,'o')
    plt.show()
    

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path

    start_time = time.time()

    train = Utils.read_train_trps_txt(path)
    print "Transaction stats"
    stats_explore(train[:,2])
    

    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])