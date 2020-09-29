# Greedy forward method
#
import networkx as nx
import random
import sys
sys.path.append("./lib")
import math
import myutilities as myutil
from Trace import Trace
import Qsap
import checkUtil
from copy import deepcopy
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
import HistoryUtil

def getUnitProb(info):
    """
    """
    return 1.0 - math.exp(-1.0*info[1][0])

def genInitialInfo(curstate,G,smodel):
    """generates initial nodes and time
    Args:
       curstate:
       G:
       smodel:
    Returns:
       inodes:   
       inittime:
    """       
    trans = {"si":Trace.S2I, "sir":Trace.S2I, "seir":Trace.S2E}
    avgspread = sum([getUnitProb(G[node1][node2][trans[smodel]]) * G[node1][node2][Trace.SPROB] for node1,node2 in G.edges()])
    avgspread /= float(G.number_of_edges())
    mylist = range(2,int(1.0/avgspread)+2)
    random.shuffle(mylist)
    inittime = int(mylist[0])
    posnodes = list(set(curstate[Trace.INFECTED]).union(curstate[Trace.RECOVERED]))
    assert len(posnodes) > 0
    random.shuffle(posnodes)
    poslen = max(1,len(posnodes)/4)
    inodes = set(posnodes[0:poslen])
    return inodes,inittime
    
  
isTest = False
def runner(inter,G,smodel,curstates,obstimes,resultfile,TESTMODE=True):
    """ runs greedy-forward baseline heuristic
    Args:
       inter:
       G:
       smodel:
       curstates:
       obstimes:
       resultfile:
       TESTMODE:
    Returns:
       history:
    """
    INFERMODE = "History"
    globals()["isTest"] = TESTMODE
    assert checkUtil.checkPreTrace(curstates,obstimes) and sorted(obstimes) == obstimes and len(obstimes) == len(curstates) and obstimes[0] == 0
    hist = {}
    if inter == "None":
       inodes,inittime = genInitialInfo(curstates[0],G,smodel)
       tracemethod = "gen{0}Trace".format(smodel.capitalize())
       method = getattr(DisTrace,tracemethod)
       (stimes,etimes,itimes,rtimes) = method(G,inodes,"static")
       trace = Trace.convertTimes2Trace(smodel,stimes,etimes,itimes,rtimes)
       hist = {0: deepcopy({Trace.SUSCEPTIBLE:set(), Trace.EXPOSED:set(), Trace.INFECTED:set(inodes),Trace.RECOVERED:set()})}
       for time in xrange(1,inittime):
           hist[time] = deepcopy(Trace.trace2Snapshot(trace,time,smodel,G))
       for tindex in xrange(0,len(curstates)-1):
           pretime,nexttime = obstimes[tindex:tindex+2]
           initnodes = curstates[tindex][Trace.INFECTED]
           if len(initnodes) == 0:
              initnodes = curstates[tindex][Trace.EXPOSED]
           assert len(initnodes) != 0    
           simlen = obstimes[tindex+1] - obstimes[tindex]
           if simlen == 1:
              continue 
           (stimes,etimes,itimes,rtimes) = method(G,initnodes,"static")
           trace = Trace.convertTimes2Trace(smodel,stimes,etimes,itimes,rtimes)
           for time in xrange(1,simlen):
               hist[time+pretime+inittime] = deepcopy(Trace.trace2Snapshot(trace,time,smodel,G))
       for obsindex in xrange(0,len(obstimes)):
           obstime = obstimes[obsindex]
           assert not hist.has_key(obstime+inittime)
           hist[obstime+inittime] = deepcopy(curstates[obsindex])
       assert len(set(hist.keys()).symmetric_difference(set(range(0,max(hist.keys())+1)))) == 0 and min(hist.keys()) == 0  
    elif inter == "bound":
       pass
    history = InputOutput.history2Times(hist,curstates,smodel) 
    if globals()["isTest"]:
       checkUtil.checkSolution(INFERMODE,inter,history,G,smodel,"notdetailed")  
    InputOutput.writeHistoryResult(history,INFERMODE,resultfile)  

    
if __name__ == "__main__":
    with gzip.open(sys.argv[1],"rb") as infile: 
        inter = cPickle.load(infile)
        G = cPickle.load(infile)
        smodel = cPickle.load(infile)
        curstates = cPickle.load(infile)
        obstimes = cPickle.load(infile)
        resultfile = cPickle.load(infile)
    runner(inter,G,smodel,curstates,obstimes,resultfile)
