#Python implementation of NetSleuth
import networkx as nx
import numpy as np
import scipy as sp
import myutilities as myutil
import scipy.sparse.linalg
from InputOutput import InputOutput
from Trace import Trace
import math
import scipy.misc

def getLN(x):
    """Ln score
    Args:
       x:
    Returns:
       score:
    """
    val = math.log(2.865064,2.0)
    last = x 
    while True:
       tscore = math.log(last,2.0)
       last = tscore
       if last > 0.0:
          val += last 
       else:
          break    
    return val

def getModelScore(S,nodenum):
    """gets model score
    Args:
       S:
       nodenum:
    Returns:
       score:
    """
    return getLN(len(S)) + math.log(scipy.misc.comb(nodenum,len(S)),2.0) 

def getP(md,fd,beta,d):
    """estimates p probability
    Args:
      md:
      fd:
      beta: unit prob
      d:
    Returns:
      score: 
    """
    pd = 1.0 - (1.0-beta)**d
    return scipy.misc.comb(fd,md)*(pd**md) * (1.0-pd)**(fd-md)


def getFperTime(Ft,md,fd,beta):
    """estimates F per time
    Args:
       Ft:
       beta:
    Returns:
       lscore:       
    """  
    lscore = 0.0
    dset = set([time])
    for d in dset:
        md = 0.0
        fd = 0.0
        lscore += getP(md,fd,beta,d) + (md*math.log(float(md)/fd)) + ((fd-md)*math.log(1.0 - float(md)/fd))
    return -1.0*lscore  
 

def getDataModelScore(S,R,G):
    """gets data given model score
    Args:
       S:
       R:
       G: 
    Returns:
       score:
    """
    avgscore = getLN(len(R.keys())) 
    for time in sorted(R.keys()):
        avgscore += getFperTime(R[time],md,fd,beta,d)
        #val = len(set(R[time]).intersection(set(S)))
        #avgscore += math.log(scipy.misc.comb(nodenum,val),2.0)
    return avgscore


def findBestRippleScore(GI,S,beta):
    """ finds the best ripple for mdl length
    Args:
       G: G only with infected nodes
       S: spreaders 
       beta: beta 
    Returns:
       attack: 
    """
    scale = int(1.0/beta)
    attack = {}
    F = {}
    curtime = 0
    while True:
       attack[curtime] = set() 
       for snode in S: 
           difset = set(GI.nodes()).difference(set(S))
           if len(difset.intersection(set(G.neighbours(snode)))) != 0:
              attack[curtime].add(snode)
    return attack

isTest = False
def runner(G,beta,inodes,resultfile,TESTMODE=False):
    """runs netsleuth
    Args:
       G:
       beta:
       inodes:
       resultfile:
       TESTMODE:
    Returns:
       spreaders: 
    """
    globals()["isTest"] = TESTMODE
    useG = nx.Graph(G)
    spreaders = set()
    lastscore = 10000000000000000.0
    while True:
       sortednodes = sorted(set(useG.nodes()).intersection(set(inodes)))
       L = nx.laplacian(useG, nodelist=sortednodes)
       eiglist,veclist = scipy.sparse.linalg.eigs(L, k=1, which='SM')
       veclist = list(veclist[:,0])
       maxindex = veclist.index(max(veclist))
       nextnode = sortednodes[maxindex]
       spreaders.add(nextnode)
       score = findBestRippleScore(useG,spreaders,beta) + getModelScore(S,useG.number_of_nodes())
       useG.remove_node(nextnode)
       dif = lastscore - score 
       if dif >= 0.0:
          break
       lastscore = score 
    history = {Trace.SUSCEPTIBLE:set(G.nodes()).difference(set(spreaders)),Trace.EXPOSED:set(),Trace.INFECTED:set(spreaders),Trace.RECOVERED:set()}, -1   
    InputOutput.writeHistoryResult(history,"Spreader",resultfile)    
    
if __name__ == "__main__":
    with gzip.open(sys.argv[1],"rb") as infile: 
        inodes = cPickle.load(infile)
        G = cPickle.load(infile)
        resultfile = cPickle.load(infile)
    runner(inodes,G,resultfile)
