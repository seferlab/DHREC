# Independent Baseline Heuristic Code
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
import HistoryUtil

def getUnitProb(info):
    """
    """
    return 1.0 - math.exp(-1.0*info[1][0])


def indepBound(node2vars,diflen,smodel):
    """runs independent baseline heuristic for the bounded case
    Args:
       node2vars:
       diflen:
       smodel:
    Returns:
       sol:
    """
    allstates = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    sol = {time: {state: set() for state in allstates} for time in xrange(diflen+1)}
    for node in node2vars.keys():
        labellist = list(node2vars[node])
        random.shuffle(labellist)
        node = int(labellist[0].replace("x","").split("?")[0])
        times = HistoryUtil.var2Times(labellist[0])
        for state in times.keys():
            sol[times[state]][state].add(node)
    return sol
                   
def indepNone(G,smodel,curstate):
    """runs independent heuristic for non bounded case
    Args:
       G:
       smodel:
       curstate:
    Returns:
       newsol:
    """
    trans = {"si":Trace.S2I, "sir":Trace.S2I, "seir":Trace.S2E}
    avgspread = sum([getUnitProb(G[node1][node2][trans[smodel]]) * G[node1][node2][Trace.SPROB] for node1,node2 in G.edges()])
    avgspread /= float(G.number_of_edges())
    mylist = range(1,int(1.0/avgspread)+1)
    random.shuffle(mylist)
    diflen = int(mylist[0])
    modcurstates = {diflen:curstate}
    node2slots = Qsap.assignStateTimeSlots(modcurstates,smodel,G)
    node2vars = Qsap.getNode2Vars(node2slots,smodel)
    sol = indepBound(node2vars,diflen,smodel)
    newsol = {time-diflen: deepcopy(sol[time]) for time in xrange(diflen)}
    return newsol

def checkIndepSolution(inter,sol):
    """checks independent solution
    Args:
       inter:
       sol:
    Returns:
       bool: true/false
    """
    if inter == "None":
       assert max(sol.keys()) < 0
    elif inter == "bound":
       assert min(sol.keys()) == 0      
    return True

isTest = False
def runner(inter,G,smodel,curstates,obstimes,resultfile,TESTMODE=True):
    """ runs independent baseline heuristics
    Args:
       inter:
       G:
       smodel:
       curstates:
       obstimes:
       resultfile:
       TESTMODE:
    Returns:
    """
    INFERMODE = "History"
    globals()["isTest"] = TESTMODE
    if inter == "None":
       assert len(curstates) == 1
       curstate = curstates[0]
       assert checkUtil.checkPreSnapshot(curstate)
       convertedhistory = indepNone(G,smodel,curstate) 
    elif inter == "bound":
       checkUtil.checkPreTrace(curstates,obstimes)
       modcurstates = {obstimes[timeindex]:deepcopy(curstates[timeindex]) for timeindex in xrange(len(obstimes))}
       node2slots = Qsap.assignStateTimeSlots(modcurstates,smodel,G)
       node2vars = Qsap.getNode2Vars(node2slots,smodel)
       diflen = max(obstimes)
       convertedhistory = indepBound(node2vars,diflen,smodel)
    print "myindo"
    print convertedhistory.keys() 
    if globals()["isTest"]:
       checkIndepSolution(inter,convertedhistory)
       checkUtil.checkSolution(INFERMODE,inter,convertedhistory,G,smodel,"notdetailed") 
    InputOutput.writeHistoryResult(convertedhistory,INFERMODE,resultfile)  

    
if __name__ == "__main__":
    with gzip.open(sys.argv[1],"rb") as infile: 
        infermode = cPickle.load(infile)
        G = cPickle.load(infile)
        smodel = cPickle.load(infile)
        curstates = cPickle.load(infile)
        obstimes = cPickle.load(infile)
        resultfile = cPickle.load(infile)
    runner(infermode,G,smodel,curstates,obstimes,resultfile)
