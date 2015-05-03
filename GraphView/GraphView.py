'''
Created on Apr 20, 2015

@author: jonathanshor
'''
import time, sys
import os.path
from optparse import OptionParser
import Utils
from Explore import Explore as ex
import numpy as np
import networkx as nx
import PerfMetrics

def get_coo(trp_raw):
# Depreciated here
# Given nX3 raw trp ndarray trp_raw, return a sparse.coo_matrix of minimal square shape
    return Utils.get_coo(trp_raw)

###############################################################################
# Update list of metrics.
###############################################################################
def append_lists(list_lists, list_vals):
    for i in range(len(list_lists)):
        list_lists[i].append(list_vals[i])

    return list_lists


def collect_degrees(coo, start_time = time.time()):
# Return 4xN ndarray of degrees, N = # of nodes
# Cols: [node, out degree, in degree, total degree]
    outs = coo.nonzero()[0]
    ins = sorted(coo.nonzero()[1])
    degs = np.union1d(outs, ins)
    degs = degs[:, np.newaxis] + [0,0,0,0]
    degs[:,1] = 0
    degs[:,2] = 0
    print "%.2f -- Collecting out-degrees" % (time.time()-start_time)
    for i in outs:
        degs[np.searchsorted(degs[:,0], i),1] += 1
    print "%.2f -- Collecting in-degrees" % (time.time()-start_time)
    for i in ins:
        degs[np.searchsorted(degs[:,0], i),2] += 1
    degs[:,3] = degs[:,1] + degs[:,2]
    print "%.2f --  avg out: %.2f, avg in: %.2f, avg tot: %.2f" % (time.time()-start_time, \
                                degs[:,1].mean(), degs[:,2].mean(), degs[:,3].mean())
    return degs

def cull_comps(G, cc_gen, giant):
# Cull G to only a giant component if giant == True, else cull the giant component
# giant comp def as a CC with more than half the vertex of entire graph
# Return culled G
    thresh = G.number_of_nodes()/2
    for cc in cc_gen:
        if len(cc) <= thresh and giant:
            # Remove small CCs from G
            G.remove_nodes_from(cc)
        if len(cc) > thresh and not giant:
            G.remove_nodes_from(cc)
    return G

def collect_comps(G, strongly, op, path):
    if strongly:
        cc_gen = nx.strongly_connected_components(G)
        ty = 'S'
    else:
        cc_gen = nx.weakly_connected_components(G)
        ty = 'W'
    if op == 1:
        ex.collect_alt_views(ex.gen_view(cc_gen), path + "%sCCsXCountView.txt" % ty, \
                             comments= "Vertex from %sCC; Count of vertex in %sCC" % (ty,ty))
    elif op == 2:
        # write raw trpl file of only vert in giant comp
        giantcc = cull_comps(G.copy(), cc_gen, True)
        fn = 'txTripletsCounts%sGiantOnly.txt' % ty
        print "Writing %s" % fn
        nx.write_weighted_edgelist(giantcc,path + fn)
        nx.write_weighted_edgelist(giantcc,'../' + fn + '.gz')
    elif op == 3:
        # write raw trpl file of only vert not in giant comp
        giantcc = cull_comps(G.copy(), cc_gen, False)
        fn = 'txTripletsCountsNo%sGiant.txt' % ty
        nx.write_weighted_edgelist(giantcc,path + fn)
    return None

def common_neighbors(G, fn, t = 0.5):
    G = G.to_undirected()
    if os.path.isfile(fn) :
        H = G.copy()
        found = nx.read_edgelist(fn, nodetype=int, data=False)
        H.add_edges_from(found.edges_iter())
        jacc_iter = nx.jaccard_coefficient(G, nx.non_edges(H))
        print "Appending to %s" % fn
        outfile = open(fn,'a',1)
        i = found.number_of_nodes()
    else:
        jacc_iter = nx.jaccard_coefficient(G)
        outfile = open(fn,'w',1)
        i = 0
    outfile.write("#vertex u; vertex v; their jaccard coef\n")
    cur = -1
    print "Starting jacc loop %s with threshold %s" % (time.strftime("%H:%M:%S"), t)
    for pair in jacc_iter:
        if pair[2] >= t:
            outfile.write("%s %s %f\n" % (pair[0],pair[1],pair[2]))
            if pair[0] != cur:
                cur = pair[0]
                i += 1
                print "%s: %s" % (i, cur)
    outfile.close()
    print "Done writing %s" % (fn)

def prune_by_degree(G, degs, k):
# Given a graph G, degree ndarray as produced by collect_degrees,
# will return G with all nodes of degree <= k removed
    prunes = degs[:,0][degs[:,3]<=k]
    G.remove_nodes_from(prunes)
    return G

def make_predG_from_jacc(undir_jaccs, dirG, testG):
# Given an undirGraph of jaccard coeffs (produced from the undirected version of dirG)
# the original digraph dirG, and the testG digraph to test against, 
# Return predicted digraph predG with arcs matching testG, each with the predicted lilekihood 
    predG = nx.DiGraph()
    for (u,v) in testG.edges_iter():
        if undir_jaccs.has_edge(u,v):
            bias = 1
            predw =undir_jaccs[u][v]['weight'] * bias 
        else:
            predw = 0
        predG.add_edge(u, v, weight=predw)
    return predG

def make_jacc_predG(dirG, testG):
# Given the original digraph dirG, and the testG digraph to test against, 
# Return predicted digraph predG with arcs matching testG,
# each with the predicted lilekihood according to Jaccard coeff
    undG = dirG.to_undirected()
    undir_jaccs = nx.Graph() 
    undir_jaccs.add_weighted_edges_from(nx.jaccard_coefficient(undG, testG.edges_iter()))
    return make_predG_from_jacc(undir_jaccs, dirG, testG)

def make_adamic_adar_index_predG(dirG, testG):
    undG = dirG.to_undirected()
    undir_AAs = nx.Graph() 
    undir_AAs.add_weighted_edges_from(nx.adamic_adar_index(undG, testG.edges_iter()))
    return make_predG_from_jacc(undir_AAs, dirG, testG)

def synced_array(arr, G):
# Given ndarray edgelist array and graph G,
# Return G as array equalling arr in columns 0 and 1
    sync = np.array(arr.copy(), dtype='float')
    for r in sync:
        r[2] = G[r[0]][r[1]]['weight']
    return sync 

def inv_shortest(pairs, G):
#     preds = np.empty_like(pairs,dtype='float') + [[0,0,0]]
    preds = np.empty((len(pairs),3),dtype='float')
    print "pairs.shope: %s, preds.shape: %s" % (pairs.shape, preds.shape)
    preds[:,0] = pairs[:,0]
    preds[:,1] = pairs[:,1]
    for r in preds:
        try:
            r[2] = 1./nx.shortest_path_length(G, int(r[0]),int(r[1]))
        except nx.NetworkXNoPath:
            r[2] = 0
    return preds

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read data from PATH', metavar='PATH')
    parser.add_option("-f", "--file", type="string", dest="fil", default='', help='file to work with')
    parser.add_option("-s", type="int", dest="strongcc", default=0, help= \
                      '1 to get SCC sizes, 2 to get trpl txt of biggest SCC, 3 to get trpl txt of all other')
    parser.add_option("-w", type="int", dest="weakcc", default=0, help= \
                      '1 to get WCC sizes, 2 to get trpl txt of biggest WCC, 3 to get trpl txt of all other')
    parser.add_option("-r", type="int", dest="prune", default=0, help= \
                      'prune all nodes with degree less than')
    parser.add_option("-d", dest="degrees", action="store_true", default=False)
    parser.add_option("-J", dest="jaccPreds", action="store_true", default=False)
    parser.add_option("-A", dest="AAPreds", action="store_true", default=False)
    parser.add_option("-P", dest="runperf", action="store_true", default=False)
    parser.add_option("-j", type="float", dest="jacc", default=-1, help='threshold for Jaccard coefs above which to record')
    (options, _args) = parser.parse_args()
    path = options.path
    print "PATH = " + path
    print "FILE = " + options.fil

    start_time = time.time()

    if options.runperf:
        list_of_metrics = [[],[],[],[],[]]
        test = Utils.read_test_trps_txt(path)
        trainx = Utils.read_train_trps_txt(path, toNX=True)

        #Inv shortest path directed
        preds = inv_shortest(test, trainx)
        (fpr, tpr, aucScore, precision, recall) = PerfMetrics.performance_metrics(test, preds[:,2])
        append_lists(list_of_metrics, [fpr, tpr, aucScore, precision, recall])
        
        #Inv shortest path reverse directed
        preds = inv_shortest(test, trainx.reverse(copy=False))
        (fpr, tpr, aucScore, precision, recall) = PerfMetrics.performance_metrics(test, preds[:,2])
        append_lists(list_of_metrics, [fpr, tpr, aucScore, precision, recall])

        #Jaccard coefficients
        preds = synced_array(test,nx.read_weighted_edgelist(path + 'jacc_preds.txt', nodetype=int))
        (fpr, tpr, aucScore, precision, recall) = PerfMetrics.performance_metrics(test, preds[:,2])
        append_lists(list_of_metrics, [fpr, tpr, aucScore, precision, recall])

        print "aucs: %s" % list_of_metrics[2]
        PerfMetrics.plotROCGraph(list_of_metrics[0],list_of_metrics[1])
        PerfMetrics.plotPRGraph(list_of_metrics[3], list_of_metrics[4])

    if options.jaccPreds:
        train = Utils.read_train_trps_txt(path, toNX=True, fn='train.txt')
        test = Utils.read_train_trps_txt(path, toNX=True, fn='test.txt')
        preds = make_jacc_predG(train, test)
        nx.write_weighted_edgelist(preds,path + 'jacc_preds.txt')

    if options.AAPreds:
        train = Utils.read_train_trps_txt(path, toNX=True, fn='train.txt')
        test = Utils.read_train_trps_txt(path, toNX=True, fn='test.txt')
        preds = make_adamic_adar_index_predG(train, test)
        nx.write_weighted_edgelist(preds,path + 'AA_preds.txt')

    if options.degrees:
        train = Utils.get_coo(Utils.read_train_trps_txt(path,fn=options.fil))
        print "%.2f -- train coo obtained." % (time.time()-start_time)
    if options.weakcc or options.strongcc or options.prune:
        trainx = Utils.read_train_trps_txt(path, toNX=True, fn=options.fil)
        print "%.2f -- nx Graph is of len: %s" % (time.time()-start_time, len(trainx))
        print nx.info(trainx)

    if options.degrees:
        ex.collect_alt_views(collect_degrees(train, start_time), path + 'DegreesView.txt', \
                             comments= "Vertex; count of out edges; count of in edges; Total adj edges")

    if options.prune:
        degs = ex.get_alt_view(path + 'DegreesView.txt')
        print "%s -- Got degs, pruning <= %s to begin" % (time.time()-start_time,options.prune)
        G = prune_by_degree(trainx, degs, options.prune)
        print "Num nodes now: %s" % G.number_of_nodes()
        nx.write_weighted_edgelist(G,path + 'PrunedBy%s.txt' % options.prune)

    if options.jacc >= 0:
        fn_base = options.fil
        trainx = nx.read_weighted_edgelist(path + fn_base, nodetype=int)
        print "%.2f -- Going for the jacc!" % (time.time() - start_time)
#         common_neighbors(trainx, path + 'JaccardCoefs', t=options.jacc)
        t=options.jacc
        common_neighbors(trainx, path + 'JaccardCoefs_t=%s' % t + fn_base, t=t)

    if options.weakcc:
        print "%.2f -- Collecting WCCs" % (time.time()-start_time)
        collect_comps(trainx, False, options.weakcc, path)
    if options.strongcc:
        print "%.2f -- Collecting SCCs" % (time.time()-start_time)
        collect_comps(trainx, True, options.strongcc, path)

    print time.time() - start_time

if __name__ == '__main__':
    main(sys.argv[1:])