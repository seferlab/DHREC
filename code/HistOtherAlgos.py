#runs other history reconstruction algos
import networkx as nx
import numpy as np
import scipy as sp
import scipy.io
import cPickle
import math
import gzip
import random
import sys
sys.path.append("./lib")
import os
from Trace import Trace
from InputOutput import InputOutput
sys.path.append("./otheralgos")
import NetSleuth
import Rumor_Centrality


def graph2Array(G):
    """converys graph to array(useful for matlab)
    Args:
       G: graph
    Returns:
       arr: array
    """
    gtype = type(G)
    maxnode = max(G.nodes())
    arr = np.zeros((maxnode+1,maxnode+1), dtype = np.int)
    for node1,node2 in G.edges():
        if gtype == nx.Graph:
           arr[node1,node2] = 1
           arr[node2,node1] = 1
        elif gtype == nx.DiGraph:   
           arr[node1,node2] = 1
    return arr

def list2str(l):
    """converts list to string for matlab
    Args:
       l: list(either single or multi dimensional)
    Returns:
       arrstr: string for matlab
    """
    lstr = " ".join([str(item) for item in l])
    return "[{0}]".format(lstr)

def arr2str(mat):
    """converts array to string for matlab
    Args:
       arr: array(either single or multi dimensional)
    Returns:
       arrstr: string for matlab
    """
    assert len(np.shape(mat)) != 1
    l1,l2 = np.shape(mat)
    matstr = ";".join([" ".join([str(mat[index1,index2]) for index2 in xrange(l2)]) for index1 in xrange(l1)])
    return "[{0}]".format(matstr[0:-1])


def runKEffectors(G,inodes,initcount):
    """ runs k effectors algorithm
    Args:
       G:
       inodes:
       initcount:
    Returns: 
       sources:    
    """
    return


def runRumorCentralityNotUsed(G,runfolder,initcount):
    """ runs Zaman's Rumor Centrality algorithm for single source node(only IC model)
    Args:
       G:
       runfolder:
       initcount: number of initial spreaders
    Returns:
       sources: initcount source nodes
    """
    codename = "rumor_centrality_all_tree"
    runfile = "temprun{0}".format(random.random())
    runpath = "{0}/{1}".format(runfolder,runfile)
    garr = graph2Array(G)
    gstr = arr2str(garr)
    code = "matlab -r \" cd('{0}'); {1}({2}, strcat('../','{3}')); cd('{4}'); quit; \" ".format(CODEDIR,codename,gstr,runpath,"..")
    os.system(code)
    result = scipy.io.loadmat(runpath)['Rlog']
    os.system("rm -rf {0}".format(runpath))
    sortedvals = sorted(result,reverse=True)[0:initcount]
    retres = [result.index(item) for item in sortedvals] 
    history = {-1: {Trace.SUSCEPTIBLE:set(G.nodes()).difference(set(spreaders)),Trace.EXPOSED:set(),Trace.INFECTED:spreaders,Trace.RECOVERED:set()}}  
    return history


def louvainCommunity(G,params):
    """fins community structure according to louvain
    Args:
       G:
       params:
    Returns:
       comdict:
    """
    CODEDIR = "louvain"
    paramstr = "-".join([part for param in params for part in param.split("/")])
    paramstr = paramstr.replace(".hist","")
    txtfile = "txt_{0}.txt".format(paramstr)
    binfile = "bin_{0}.bin".format(paramstr)
    treefile = "tree_{0}.tree".format(paramstr)
    if G.number_of_nodes() == 1:
       comdict = {G.nodes()[0]: 1} 
       return comdict  
    node2index,index2node = {},{}
    allnodes = G.nodes()
    for index in xrange(G.number_of_nodes()):
        node2index[allnodes[index]] = index
        index2node[index] = allnodes[index]
    with open(txtfile,"w") as outfile:
        for node1,node2 in G.edges():
            outfile.write("{0}\t{1}\n".format(node2index[node1],node2index[node2]))    
    code =  "{0}/convert -i {1} -o {2}".format(CODEDIR,txtfile,binfile)
    os.system(code)
    code =  "{0}/community {1} -l -1 > {2}".format(CODEDIR,binfile,treefile)
    os.system(code)
    comdict = {}
    with open(treefile,"r") as infile:
        for line in infile:
            line = line.rstrip()
            node,com = [int(part) for part in line.split(" ")]
            if comdict.has_key(index2node[node]):
               break 
            comdict[index2node[node]] = com            
    assert len(set(comdict.keys()).symmetric_difference(G.nodes())) == 0
    for filename in [txtfile,binfile,treefile]:
        delcode = "rm -rf {0}".format(filename)
        os.system(delcode)
    return comdict
  

def runOther(graphfile,algo,algoparams,snapshotfile,resultfile,smodel,runfolder,inter):    
    """ runs Other Algos
    Args:  
       graphfile:
       algo:
       algoparams:
       snapshotfile:
       resultfile:
       smodel:
       runfolder:
       inter:
    """
    G = InputOutput.readGraphAndParams(graphfile,outformat="pickle")     
    curinfo,obstimes = InputOutput.readSnapshot(snapshotfile,"dis")
    curinfo = curinfo[0]
    assert len(obstimes) == 1 and algo in ["NetSleuth","RumorCentrality","KEffectors"]
    methodname = "run{0}".format(algo)
    if algo == "NetSleuth":
       pureG = nx.DiGraph()
       avgs,avgpar = 0.0, 0.0
       for node1,node2 in G.edges():
           pureG.add_edge(node1,node2)
           avgs += G[node1][node2][Trace.SPROB]
           if smodel in ["si","sir"]:
              avgpar += G[node1][node2][Trace.S2I][1][0]
           elif smodel == "seir":
              avgpar += G[node1][node2][Trace.S2E][1][0]     
       avgs /= float(G.number_of_edges())
       avgpar /= float(G.number_of_edges())
       beta = avgs * (1.0 - math.e**(-1.0*avgpar))
       inodes = list(curinfo[Trace.INFECTED])
       NetSleuth.runner(pureG,beta,inodes,resultfile,True)
    elif algo == "RumorCentrality":
       if len(curinfo[Trace.INFECTED]) == 0:
          return  
       subG = nx.Graph(G.subgraph(curinfo[Trace.INFECTED]))
       allnodes = subG.nodes()
       while nx.number_connected_components(subG) != 1:
          random.shuffle(allnodes)
          node1,node2 = allnodes[0:2]
          if not subG.has_edge(node1,node2):
             subG.add_edge(node1,node2)
       assert nx.number_connected_components(subG) == 1
       comdict = louvainCommunity(subG,params=[resultfile])
       initcount = len(set(comdict.values()))
       if smodel == "si":
          inodes = list(curinfo[Trace.INFECTED])
       elif smodel in ["sir","seir"]:
          inodes = list(curinfo[Trace.INFECTED]) + list(curinfo[Trace.RECOVERED])
       if len(inodes) < initcount:
          initcount = len(inodes)
       Rumor_Centrality.runner(G,initcount,inodes,resultfile,True) 
    elif algo == "KEffectors":
       subG = nx.Graph(G.subgraph(curinfo[Trace.INFECTED]))
       comdict = louvainCommunity(G)
       initcount = len(set(com for com in comdict.values()))
       print "Not implemented yet!!"
       exit(1) 
       inodes = list(curinfo[Trace.INFECTED])
       history = globals()[methodname](G,inodes,runfolder)
       InputOutput.writeHistoryResult(history,"Spreader",resultfile)           
    

MATLABPATH = "/usr/local/bin/matlab"
CODEDIR = "otheralgos"

def main():
    with gzip.open(sys.argv[1],"rb") as file:
        graphfile = cPickle.load(file)
        algo = cPickle.load(file)
        algoparams = cPickle.load(file)
        snapshotfile = cPickle.load(file)
        resultfile = cPickle.load(file)
        smodel = cPickle.load(file)
        trunfolder = cPickle.load(file)
        inter =  cPickle.load(file)
    runOther(graphfile,algo,algoparams,snapshotfile,resultfile,smodel,trunfolder,inter)

if __name__ == "__main__":
    main()
