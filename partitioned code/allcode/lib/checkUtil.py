#Methods for testing
import random
import sys
sys.path.append("./lib")
import networkx as nx
import myutilities as myutil
import os
import HistoryUtil
import math
import itertools
from Trace import Trace


def checkPreSnapshot(curstate):
    """checks snapshot data before running method
    Args:
       curstate: given snapshot
    Returns:
       bool: true/false
    """
    impnodes = set(curstate[Trace.INFECTED]).union(curstate[Trace.RECOVERED]).union(curstate[Trace.EXPOSED])
    if len(impnodes) == 0:
       print "No data is given for inference"
       return False
    return True

def checkPreTrace(curstates,obstimes):
    """checks trace data before algo execution
    Args:
       curstates:
       obstimes:
    Returns:
       bool: true/false
    """
    if len(obstimes) == 0:
       print "No data given for inference"
       return False
    assert obstimes == sorted(obstimes)
    for curstate in curstates:
        if checkPreSnapshot(curstate):
           return True 
    print "Given Diffusion Data is not good!!"    
    return False

def checkTransitionCount(sol,smodel):
    """checks state transition count in solution format
    Args:
      sol:
      smodel:
    Returns:
      bool: true/false 
    """
    allstates = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    nodetimes = {state: [] for state in allstates}
    for time in sol.keys():
        for state in sol[time].keys():
            nodetimes[state].extend(list(sol[time][state]))
    for state in nodetimes.keys():
        assert len(nodetimes[state]) == len(set(nodetimes[state]))        


def his2StateTimes(hist,state):
    """extract times from hist
    Args:
       hist:
       state:
    Returns:
       times:
    """
    times = {}
    for time in hist.keys():
        if not hist.has_key(state):
           continue 
        for node in hist[state]:
            assert not times.has_key(node)
            times[node] = time        
    return times
    
        
def checkCover(hist,G,smodel):
    """checks cover constraints
    Args:
       hist: history of transition times
       G:
       smodel:
    Returns:
       truth: true/false
    """
    etimes = his2StateTimes(hist,Trace.EXPOSED)
    itimes = his2StateTimes(hist,Trace.INFECTED)
    rtimes = his2StateTimes(hist,Trace.RECOVERED)
    mintime = min(etimes.values()+itimes.values()+rtimes.values())
    if smodel == "seir":
       rectimes, sendtimes = etimes,itimes 
    elif smodel in ["si","sir"]:
       rectimes, sendtimes = itimes,itimes
    startnodes = set(node for node in sendtimes.keys() if sendtimes[node] == mintime)     
    for recnode in rectimes.keys():
        if recnode in startnodes:
           continue 
        valid = False
        for sendnode in G.predecessors(recnode):
            if G.has_edge(sendnode,recnode) and sendtimes.has_key(sendnode) and sendtimes[sendnode] < rectimes[recnode] and (not rtimes.has_key(sendnode) or (rtimes.has_key(sendnode) and rtimes[sendnode] >= rectimes[recnode])):
               valid = True
               break
        if not valid:
           return False
    return True


def checkCoverBetween(prestate,curstate,smodel,G):
    """given all state info between two time points, checks covering constraints
    Args:
       prestate: state at j-1
       curstate: state at j
       smodel:
       G:
    Returns:
       bool: true/false
    """
    assert len([node for state in prestate.keys() for node in prestate[state]]) == len([node for state in prestate.keys() for node in prestate[state]])
    if smodel == "seir":
       giver,taker = Trace.INFECTED, Trace.EXPOSED
    elif smodel in ["si","sir"]:    
       giver,taker = Trace.INFECTED, Trace.INFECTED 
    newtakenodes = curstate[taker].difference(set(prestate[taker]))  
    print prestate[giver]
    print curstate[taker]  
    for takenode in newtakenodes:
        print takenode,len(set(G.predecessors(takenode)).intersection(prestate[giver]))
        assert len(set(G.predecessors(takenode)).intersection(prestate[giver])) > 0     
    return True    

        
def checkDiffusionOrder(transtimes,smodel,G):
    """check order of diffusion progression
    Args:
       transtimes: history in terms of exact transition times (time, state, nodes)
       smodel:
       G:
    Returns:
       bool: true/false
    """
    mintime = min(transtimes.keys())
    etimes,itimes,rtimes = HistoryUtil.getTimes(transtimes,Trace.EXPOSED), HistoryUtil.getTimes(transtimes,Trace.INFECTED), HistoryUtil.getTimes(transtimes,Trace.RECOVERED)
    if smodel == "seir":
       sender,receiver = itimes,etimes
    elif smodel in ["si","sir"]: 
       sender,receiver = itimes,itimes    
    for node in receiver.keys():
        if receiver[node] != mintime:
           flag = False 
           for neighnode in G.predecessors(node):
               if sender.has_key(neighnode) and sender[neighnode] < receiver[node]:
                  if smodel in ["sir","seir"] and rtimes.has_key(neighnode) and rtimes[neighnode] < receiver[node]:
                     continue  
                  flag = True
                  break 
           assert flag
    return True
    
        
def checkSolution(infermode,inter,history,G,smodel,check="detailed"):
    """tests whether solution is valid or not(solution format)
    Args:
       infermode:
       inter:
       history:
       G:
       smodel:
       check: sol is valid or not
    Returns:
       bool: true/false
    """   
    if smodel == "seir":
       giver,taker = Trace.INFECTED, Trace.EXPOSED
    elif smodel in ["si","sir"]:    
       giver,taker = Trace.INFECTED, Trace.INFECTED   
    allstates = [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    if infermode == "History": 
       for time in history.keys():
           assert len(set(allstates).difference(set(history[time].keys()))) == 0 
       mintime = min(history.keys())     
       if inter != "bound" or check == "detailed":
          assert len(history[mintime][Trace.INFECTED]) != 0
       for state in set(history[mintime].keys()).difference(set([Trace.SUSCEPTIBLE,Trace.INFECTED])):
           print state,len(history[mintime][state]) 
           assert len(history[mintime][state]) == 0
       soltimes = history.keys()
       checkTransitionCount(history,smodel)
       if check == "detailed":
          checkDiffusionOrder(history,smodel,G)
    elif infermode == "Spreader":
       if inter in ["Bound","None"]: 
          assert len(history[0][Trace.EXPOSED]) == 0 and len(history[0][Trace.RECOVERED]) == 0 and len(history[0][Trace.INFECTED]) != 0  
    return True


def basicCheck(history,G):
    """history in state format(basic check)
    Args:
       history:
       G:
    Returns:
       bool: true/false
    """
    for time in history.keys():
        print time,len(history[time][Trace.SUSCEPTIBLE]),len(history[time][Trace.EXPOSED]), len(history[time][Trace.INFECTED]),len(history[time][Trace.RECOVERED])
        assert basicStateCheck(history[time],G)
    return True

def basicStateCheck(curstate,G):
    """basic check for state
    Args:
       curstate:
       G:
    Returns:
       bool: true/false
    """
    assert sum([len(curstate[state]) for state in curstate.keys()]) == G.number_of_nodes()
    for state1,state2 in list(itertools.product(curstate.keys(),curstate.keys())):
        if state1 != state2:
           assert len(curstate[state1].intersection(curstate[state2])) == 0
    return True


def isSnapshotValid(prestate,curstate,G,smodel,check="detailed"):
    """check whether returned snapshot is valid or not(not solution format)
    Args:
       prestate:
       curstate:
       G:
       smodel:
       check: whether check is detailed or not
    Returns:
       bool: true/false
    """
    assert basicStateCheck(curstate,G) and basicStateCheck(prestate,G)
    assert len(set(curstate[Trace.SUSCEPTIBLE]).difference(set(prestate[Trace.SUSCEPTIBLE]))) == 0
    assert len(set(prestate[Trace.RECOVERED]).difference(set(curstate[Trace.RECOVERED]))) == 0 
    if check == "detailed":
       checkCoverBetween(prestate,curstate,smodel,G)
    return True


def checkPreConversion(sol,G,smodel):
    """checks solution before converting to solution format
    Args: 
       sol:
       G:
       smodel:
    Returns:
       bool: true/false
    """
    soltimes = sorted(sol.keys())
    mintime = min(soltimes)
    for tindex in xrange(1,len(soltimes)):
        pretime,curtime = soltimes[tindex-1:tindex+1]
        prestate,curstate = sol[pretime], sol[curtime]
        isSnapshotValid(prestate,curstate,G,smodel,"notdetailed")
    assert len(sol[mintime][Trace.INFECTED]) != 0 and len(sol[mintime][Trace.EXPOSED]) == 0 and len(sol[mintime][Trace.RECOVERED]) == 0      
    return True
  

def checkStateOrder(hist,smodel):
    """checks ordering among rounded vars
    Args:
       hist:
       smodel:
    Returns:
       bool: true/false
    """
    times = {state: HistoryUtil.getTimes(hist,state) for state in [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]}
    allnodes = set([node for time in hist.keys() for node in hist[time].keys()])
    mintime = min(hist.keys())
    startnodes = set(node for node in itimes.keys() if itimes[node] == mintime)
    if smodel == "si":
       statepairs = []
    elif smodel == "sir":
       statepairs = [(Trace.INFECTED, Trace.RECOVERED)]
    elif smodel == "seir":
       statepairs = [(Trace.EXPOSED, Trace.INFECTED),(Trace.INFECTED, Trace.RECOVERED)]
    for node in allnodes:
        for state1,state2 in statepairs:
            if not (state1 == Trace.EXPOSED and state2 == Trace.INFECTED and smodel == "seir") and not times[state1].has_key(node) and times[state2].has_key(node):
               return False 
            if times[state1].has_key(node) and times[state2].has_key(node) and times[state1][node] >= times[state2][node]:
               return False 
    return True


def checkStateExistence(hist,seenstate,smodel):
    """checks state existence
    Args:
       hist:
       seenstate:
       smodel:
    Returns:
       bool: true/false
    """
    etimes,itimes,rtimes = HistoryUtil.getTimes(hist,Trace.EXPOSED), HistoryUtil.getTimes(hist,Trace.INFECTED), HistoryUtil.getTimes(hist,Trace.RECOVERED)
    mintime = min(hist.keys())
    startnodes = set(node for node in itimes.keys() if itimes[node] == mintime) 
    for node in seenstate[Trace.EXPOSED]:
        if not etimes.has_key(node):
           return False 
    for node in seenstate[Trace.INFECTED]:
        if not itimes.has_key(node) or (smodel == "seir" and node not in startnodes and not etimes.has_key(node)):
           return False 
    for node in seenstate[Trace.RECOVERED]:
        if not rtimes.has_key(node) or not itimes.has_key(node) or (smodel == "seir" and node not in startnodes and not etimes.has_key(node)):
           return False 
    return True        
            
def isHistoryValid(hist,smodel,G,seenstate,noisy):
    """checks wheter history is valid or not
    Args:
       hist:
       smodel:
       G:
       seenstate:
       noisy:
    Returns:
       bool: true/false
    """
    params = [etimes,itimes,rtimes] 
    return checkStateExistence(hist,seenstate,smodel) and checkCover(hist,G,smodel,params) and checkStateOrder(hist,smodel,params)
