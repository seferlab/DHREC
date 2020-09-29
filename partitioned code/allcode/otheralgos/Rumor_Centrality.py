#Python implementation of Rumor Centrality
import networkx as nx
import myutilities as myutil
from InputOutput import InputOutput
from Trace import Trace
import math
import operator
import random

def estimateScores(tree,source,t,plog):
    """estimates rumor score given t p
    Args:
       tree:
       source:
       t:
       plog:
    Returns:
       logscore:
    """
    N = tree.number_of_nodes()
    if N > 10:  #use Stirling approx for log N!
       z = N*(math.log(N)-1) + (0.5*math.log(2.0*math.pi*N))
    else:
       z = math.log(math.factorial(N))
    logscore = {node:None for node in tree.nodes()}
    sortednodes = nx.dfs_preorder_nodes(tree,source)
    for node in sortednodes:
        if node == source:
           logscore[node] = z - math.log(N) - sum([plog[child] for child in tree.successors(node)]) 
        else:
           parent = tree.predecessors(node)[0] 
           assert logscore[parent] != None
           logscore[node] = logscore[parent] + math.log(t[node]) - math.log(N-t[node])
    for node in tree.nodes():
        assert logscore[node] != None
    return logscore 


def estimateLogCount(tree,source):
    """estimates log counts t,p
    Args: 
       tree:
       source:
    Returns:
       t:
       plog:
    """
    t = {node: None for node in tree.nodes() if node != source}
    plog = {node: None for node in tree.nodes() if node != source}
    sortednodes = nx.dfs_postorder_nodes(tree,source)
    for node in sortednodes:
        if len(tree.successors(node)) == 0:
           t[node] = 1
           plog[node] = 0
        elif node != source:
           t[node] = sum([t[child] for child in tree.successors(node)]) + 1
           plog[node] = math.log(t[node]) + sum([plog[child] for child in tree.successors(node)])
    assert not t.has_key(source) and not plog.has_key(source)       
    for val in t.values():
        assert val != None
    for val in plog.values():
        assert val != None
    return t,plog

def findBfsTree(G):
    """returns bfs tree of G
    Args:
       G:
    Returns:
       tree:
       source:   
    """  
    bfsfound,findsource = False, None
    for source in G.nodes():
        flag = True
        for target in G.nodes():
            if source == target:
               continue 
            if not nx.has_path(G, source, target):
              flag = False
              break
        if flag:    
           bfsfound = True 
           tree = nx.bfs_tree(G,source)
           findsource = source
           break
    if not bfsfound:
       print "ERROR: No bfs tree defined over G"
       exit(1) 
    return tree,findsource


isTest = False
def runner(G,initcount,inodes,resultfile,TESTMODE=False):
    """runs rumor centrality
    Args:
       G: 
       initcount: number of initial spreaders
       inodes:
       resultfile:
       TESTMODE:
    Returns:
       spreaders: 
    """
    globals()["isTest"] = TESTMODE
    G = nx.Graph(G)
    subG = G.subgraph(inodes)
    allnodes = subG.nodes()
    while nx.number_connected_components(subG) != 1:
       random.shuffle(allnodes)
       node1,node2 = allnodes[0:2]
       if not subG.has_edge(node1,node2):
          subG.add_edge(node1,node2)
    assert nx.number_connected_components(subG) == 1		
    tree,source = findBfsTree(subG)
    t,plog = estimateLogCount(tree,source)
    logrumor = estimateScores(tree,source,t,plog)  
    spreaders = [node for node,val in sorted(logrumor.iteritems(), key=operator.itemgetter(1), reverse=True)[0:initcount]]
    history = {Trace.SUSCEPTIBLE:set(G.nodes()).difference(set(spreaders)),Trace.EXPOSED:set(),Trace.INFECTED:set(spreaders),Trace.RECOVERED:set()}, -1   
    InputOutput.writeHistoryResult(history,"Spreader",resultfile)  
       
if __name__ == "__main__":
    with gzip.open(sys.argv[1],"rb") as infile:
        G = cPickle.load(infile)
        initcout = cPickle.load(infile)
        inodes = cPickle.load(infile)
        resultfile = cPickle.load(infile)
    runner(G,initcount,inodes,resultfile)
