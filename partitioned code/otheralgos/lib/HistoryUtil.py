"""Diffusion History Reconstruction Utilities related to History Inference
"""
import networkx as nx
import os
import sys
import math
import random
import time
import numpy as np
import myutilities as myutil
from Trace import Trace


def getStates(smodel):
    """returns states
    Args:
       smodel:
    Returns:
       states:   
    """
    states = []
    if smodel in ["si","sis"]:
       states = [Trace.SUSCEPTIBLE,Trace.INFECTED]
    elif smodel == "sir":
       states = [Trace.SUSCEPTIBLE,Trace.INFECTED,Trace.RECOVERED]
    elif smodel == "seir":
       states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    return states

    
def removeZeros(retvalues):
    """removes zeros
    Args:
       retvalues:
    Returns:
    """
    retkeys = retvalues.keys()    
    for key in retkeys:
        if retvalues[key] in [0.0,-0.0]:
           del retvalues[key] 
    return retvalues


def runHistoryCode(consstr,objstr,boundstr,runfolder,varstarter,filestr,lastmethod=None):
    """runs cplex code from history constraints
    Args:
        consstr: constraint string
        objstr: objective function string
        boundstr: boundary string
        runfolder: run folder
        varstarter: variable start symbol ["x"]
        filestr: filename str
        lastmethod: method to run before returning
    Returns:
        edge2val : edge value mapping
    """
    PNUM = 1
    outlppath = "{0}/{1}.lp".format(runfolder,filestr)
    with open(outlppath,"w") as file:
       file.write(objstr+"\n")
       file.write(consstr+"\n")
       file.write(boundstr+"\n")
       file.write("End\n")
    cplexoutpath = "{0}/{1}.lpout".format(runfolder,filestr)
    cplexscriptpath = "{0}/{1}.script".format(runfolder,filestr)
    with open(cplexscriptpath,"w") as file:
       file.write("read {0}\n".format(outlppath))
       file.write("set threads {0}\n".format(PNUM))
       file.write("optimize\n") 
       file.write("display solution objective\n")
       file.write("display solution variables -\n")  
    t1 = time.time()
    code = "cplex < {0} > {1}".format(cplexscriptpath,cplexoutpath)
    os.system(code)
    t2 = time.time()    
    runtime = t2-t1
    from InputOutput import InputOutput
    retvalues,objval = InputOutput.readCplexOut(cplexoutpath,specific = varstarter)
    assert checkValidCplexRun(cplexoutpath)
    for cplexfile in myutil.listfiles(runfolder):
        os.system("rm -rf {0}/{1}".format(runfolder,cplexfile))
    if lastmethod != None:
       retvalues = globals()[lastmethod](retvalues) 
    return retvalues,objval,runtime



def checkAffect(sendtimes,rectimes,smodel,sendstate,recstate,info="normal"):
    """check whether var1 affects var2(whether var1 is before var2) does not consider recovered info!! 
    Args:
       sendtimes:
       rectimes:
       smodel:
       sendstate:
       recstate:
       degree:
    Returns:
       bool: true/false
    """
    assert info in ["weak","normal"]
    if info == "weak":
       if sendtimes.has_key(sendstate) and rectimes.has_key(recstate): 
          return sendtimes[sendstate] < rectimes[recstate]
       return False
    elif info == "normal":
       if sendtimes.has_key(sendstate) and rectimes.has_key(recstate): 
          return sendtimes[sendstate] < rectimes[recstate] and (not sendtimes.has_key(Trace.RECOVERED) or (sendtimes.has_key(Trace.RECOVERED) and sendtimes[Trace.RECOVERED] >= rectimes[recstate]))
       return False 


def var2Times(varname):
    """extract times from varname
    Args:
       varname:
    Returns:
       times:
    """            
    splitted = varname.replace("x","").split("?")[1:]
    times = {splitted[2*index+0]: int(splitted[2*index+1]) for index in xrange(len(splitted)/2)}
    return times


def buildSnapshotFromTransTimes(transtimes,smodel,G):
    """builds snapshot type from transtimes
    Args:
       transtimes: 
       smodel:
       G:
    Returns
       snaphist:
    """
    states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    allnodes = G.nodes()
    snaphist = {time: {state:set() for state in states} for time in transtimes.keys()}
    mintime, maxtime = min(snaphist.keys()), max(snaphist.keys())
    INFTIME = maxtime+1
    if smodel == "si":
       itimes = getTimes(transtimes,Trace.INFECTED)
       for node in allnodes:
           itime = INFTIME
           if itimes.has_key(node):
              itime = itimes[node]
           for tindex in xrange(mintime,itime):
               snaphist[tindex][Trace.SUSCEPTIBLE].add(node)
           for tindex in xrange(itime,maxtime+1):
               snaphist[tindex][Trace.INFECTED].add(node)
    elif smodel == "sir":
       itimes = getTimes(transtimes,Trace.INFECTED)
       rtimes = getTimes(transtimes,Trace.RECOVERED)
       for node in allnodes:
           itime,rtime = INFTIME, INFTIME
           if itimes.has_key(node):
              itime = itimes[node]   
           if rtimes.has_key(node):
              rtime = rtimes[node]
           for tindex in xrange(mintime,itime):
               snaphist[tindex][Trace.SUSCEPTIBLE].add(node)
           for tindex in xrange(itime,rtime):
               snaphist[tindex][Trace.INFECTED].add(node)
           for tindex in xrange(rtime,maxtime+1):
              snaphist[tindex][Trace.RECOVERED].add(node)   
    elif smodel == "seir":
       etimes = getTimes(transtimes,Trace.EXPOSED)
       itimes = getTimes(transtimes,Trace.INFECTED)
       rtimes = getTimes(transtimes,Trace.RECOVERED)
       for node in allnodes:
           etime,itime,rtime = INFTIME,INFTIME, INFTIME
           if etimes.has_key(node):
              etime = etimes[node]   
           if itimes.has_key(node):
              itime = itimes[node]   
           if rtimes.has_key(node):
              rtime = rtimes[node]
           if etime != INFTIME:
              for tindex in xrange(mintime,etime):
                  snaphist[tindex][Trace.SUSCEPTIBLE].add(node)
              for tindex in xrange(etime,itime):
                  snaphist[tindex][Trace.EXPOSED].add(node)
           else:  
              for tindex in xrange(mintime,itime):
                  snaphist[tindex][Trace.SUSCEPTIBLE].add(node)
           for tindex in xrange(itime,rtime):
               snaphist[tindex][Trace.INFECTED].add(node)
           for tindex in xrange(rtime,maxtime+1):
               snaphist[tindex][Trace.RECOVERED].add(node)    
    return snaphist


def assignStateTransitions(curstate,prestate,smodel):
    """assigns state transitions
    Args:
       curstate:
       prestate:
    Returns: 
       transnodes:   
    """
    tstates = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    statepairs = []
    for index in xrange(len(tstates)):
        statepairs.append((tstates[index],tstates[index]))
        if index != len(tstates) - 1:
           statepairs.append((tstates[index],tstates[index+1]))
    transnodes = {}    
    for state1,state2 in statepairs:
        transnodes[(state1,state2)] = set(prestate[state1]).intersection(set(curstate[state2]))
    return transnodes


def checkValidCplexRun(cplexoutpath):
    """checks whether cplex run produce any errpr
    Args: 
       cplexoutpath:
    Returns:
       bool: false/true
    """
    with open(cplexoutpath,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if line.find("CPXERR") != -1 or line.find("ERROR") != -1 or line.find("Error") != -1 or line.find("No problem exists") != -1:
               return False 
    return True


def FixSameParts(sol,G):
    """only delete equal parts when inter == None
    Args:
       sol: 
       G:
    Returns:
       sol: 
    """
    mintime = min(sol.keys())
    sametimes = set([mintime])  
    for time in sorted(sol.keys()):
        flag = True
        for state in sol[time].keys():
            if len(sol[time][state].difference(sol[mintime][state])) != 0 or len(sol[mintime][state].difference(sol[time][state])) != 0:
               flag = False
               break
        if not flag:
           break
        else:
           sametimes.add(time)
    maxsametime = max(sametimes)       
    for time in sametimes.difference(set([maxsametime])):
        del sol[time]
    if len(sol.keys()) == 1:
       return sol        
    allstimes = set()  
    for time in sorted(sol.keys()):
        if len(sol[time][Trace.SUSCEPTIBLE]) == G.number_of_nodes():
           allstimes.add(time)
        else:
           break  
    for time in allstimes:
        del sol[time]    
    return sol    
        
def getTimes(hist,state):
    """ get state times from hist type(exact state transition times)
    Args:
       hist:
       state:
    Returns:
       times:
    """
    times = {}
    for time in hist.keys():
        if hist[time].has_key(state):
           for node in hist[time][state]:
               times[node] = time 
    return times


def lastFixSolution(sol,G,smodel,inter,infermode):   
    """makes last fixes to solution
    Args:
       sol: current solution
       G: Graph
       smodel: spreading model
       inter: whether bound is given or not
       infermode:
    Returns:
       sol:  
    """  
    if inter == "None":
       return lastFixNone(sol,G,smodel,infermode)
    elif inter == "bound":
       return lastFixBound(sol,G,smodel,infermode)

def lastFixBound(sol,G,smodel,infermode): 
    """makes last solution fixes in bound case
    Args:
       sol: current solution
       G: Graph
       smodel: spreading model 
       infermode: 
    Returns:
       sol: 
    """
    if infermode == "History":
       return sol
    elif infermode == "Spreader":  
       mintime = None
       for time in sorted(sol.keys()):
           if len(sol[time][Trace.INFECTED]) != 0:
              mintime = time
              break 
       assert mintime != None
       if smodel == "sir":
          if len(sol[mintime][Trace.RECOVERED]) != 0:
             rsnodes = set(node for node in sol[mintime][Trace.RECOVERED] if getUnitProb(G.node[node][Trace.I2R]) >= random.random()) 
             rinodes = sol[mintime][Trace.RECOVERED].difference(rsnodes)
             if len(rinodes) == 0 and len(sol[mintime][Trace.INFECTED]) == 0:
                rinodes = set([random.choice(list(rsnodes))])
                rsnodes = rsnodes.difference(rinodes) 
             sol[mintime][Trace.INFECTED] |= rinodes
             sol[mintime][Trace.SUSCEPTIBLE] |= rsnodes
             sol[mintime][Trace.RECOVERED] = set() 
       elif smodel == "seir":
          if len(sol[mintime][Trace.EXPOSED]) != 0:
             rsnodes = sol[mintime][Trace.EXPOSED]
             sol[mintime][Trace.SUSCEPTIBLE] |= rsnodes #?
             sol[mintime][Trace.EXPOSED] -= rsnodes
          elif len(sol[mintime][Trace.RECOVERED]) != 0:
             rsnodes = set(node for node in sol[mintime][Trace.RECOVERED] if getUnitProb(G.node[node][Trace.I2R]) >= random.random()) 
             rinodes = sol[mintime][Trace.RECOVERED].difference(set(rsnodes))
             if len(rinodes) == 0 and len(sol[mintime][Trace.INFECTED]) == 0:
                rinodes = set([random.choice(list(rsnodes))])
                rsnodes = rsnodes.difference(rinodes) 
             sol[mintime][Trace.INFECTED] |= rinodes
             sol[mintime][Trace.SUSCEPTIBLE] |= rsnodes
             sol[mintime][Trace.RECOVERED] = set()  
    return sol   

def getUnitProb(info):
    """
    """
    return 1.0 - math.exp(-1.0*info[1][0])
    
def lastFixNone(sol,G,smodel,infermode): 
    """makes last solution fixes in NONE case
    Args:
       sol: current solution
       G: Graph
       smodel: spreading model 
       infermode:
    Returns:
       sol: 
    """
    sol = FixSameParts(sol,G)
    mintime = min(sol.keys())
    print "solkeys:"
    print sol.keys()
    if smodel == "sir":
       if len(sol[mintime][Trace.RECOVERED]) != 0:
          rsnodes = set(node for node in sol[mintime][Trace.RECOVERED] if getUnitProb(G.node[node][Trace.I2R]) >= random.random()) 
          rinodes = sol[mintime][Trace.RECOVERED].difference(rsnodes)
          if len(rinodes) == 0 and len(sol[mintime][Trace.INFECTED]) == 0:
             rinodes = set([random.choice(list(rsnodes))])
             rsnodes = rsnodes.difference(rinodes)      
          sol[mintime][Trace.INFECTED] |= set(rinodes)
          sol[mintime][Trace.SUSCEPTIBLE] |= set(rsnodes)
          sol[mintime][Trace.RECOVERED] = set()
          for node in rinodes:
              rtime = mintime+1    
              for time in xrange(min(sol.keys())+1,rtime):
                  sol[time][Trace.RECOVERED].remove(node)
          for node in rsnodes:
              rtime = mintime+2
              itime = mintime+1    
              for time in xrange(itime,rtime):
                  if sol.has_key(time):
                     sol[time][Trace.INFECTED].add(node)  
              for time in xrange(min(sol.keys())+1,rtime):  
                  if sol.has_key(time):    
                     sol[time][Trace.RECOVERED].remove(node)
    elif smodel == "seir":
       if len(sol[mintime][Trace.EXPOSED]) != 0:
          esnodes = sol[mintime][Trace.EXPOSED]
          sol[mintime][Trace.INFECTED] |= esnodes 
          sol[mintime][Trace.EXPOSED] -= esnodes
       if len(sol[mintime][Trace.RECOVERED]) != 0:
          rsnodes = set(node for node in sol[mintime][Trace.RECOVERED] if getUnitProb(G.node[node][Trace.I2R]) >= random.random()) 
          rinodes = sol[mintime][Trace.RECOVERED].difference(set(rsnodes))
          if len(rinodes) == 0 and len(sol[mintime][Trace.INFECTED]) == 0:
             rinodes = set([random.choice(list(rsnodes))])
             rsnodes = rsnodes.difference(rinodes) 
          sol[mintime][Trace.INFECTED] |= rinodes
          sol[mintime][Trace.SUSCEPTIBLE] |= rsnodes
          sol[mintime][Trace.RECOVERED] = set() 
          for node in rinodes:
              rtime = mintime+1    
              for time in xrange(min(sol.keys())+1,rtime):
                  sol[time][Trace.RECOVERED].remove(node)
          for node in rsnodes:
              etime = mintime+1
              itime = mintime+2
              rtime = mintime+3    
              for time in xrange(etime,itime):
                  if sol.has_key(time):
                     sol[time][Trace.EXPOSED].add(node)  
              for time in xrange(itime,rtime):
                  if sol.has_key(time):
                     sol[time][Trace.INFECTED].add(node)      
              for time in xrange(min(sol.keys())+1,rtime):
                  if sol.has_key(time):      
                     sol[time][Trace.RECOVERED].remove(node)
       print "after modi:"              
       for time in sorted(sol.keys()):
           print time,len(sol[time][Trace.SUSCEPTIBLE]),len(sol[time][Trace.EXPOSED]),len(sol[time][Trace.INFECTED]),len(sol[time][Trace.RECOVERED])
      
    return sol


def isPerfectStopCond(history,G,smodel,typemethod):
    """should we stop inferring for perfect case
       Args:
         history: all history inferred so far
         G: Graph
         smodel: spreading model
         typemethod: type method
       Returns:
         flag: boolean: should algorithm stop or not True->stop, False->not stop
    """
    assert typemethod in ["single","double","Pcdsvc","Pcvc","MinCut"]
    if typemethod in ["single","double","Pcdsvc","Pcvc","MinCut"]:
       STEPCOUNT = 5
       if len(history) < STEPCOUNT:
          return False
       mintime = min(history.keys()) 
       #if len(history[mintime][Trace.SUSCEPTIBLE]) == G.number_of_nodes():
       #   return True
       for index in xrange(STEPCOUNT-1):
           for state in history[mintime + index].keys():
               if len(set(history[mintime + index + 1][state]).difference(set(history[mintime + index][state]))) != 0 or len(set(history[mintime + index][state]).difference(set(history[mintime + index + 1][state]))) != 0:
                  return False
    elif typemethod == "frac":
       FRACEPSILON = 1.0 
       mysum = 0.0
       STEPCOUNT = 2
       if len(history) < STEPCOUNT:
          return False
       mintime = min(history.keys())
       for statekey in history[mintime].keys():
           dif1 = set(history[mintime][statekey].keys()).difference(set(history[mintime + 1][statekey].keys()))
           inter = set(history[mintime][statekey].keys()).intersection(set(history[mintime + 1][statekey].keys()))
           dif2 = set(history[mintime + 1][statekey].keys()).difference(set(history[mintime][statekey].keys()))
           mysum += sum([0.0] + [ abs(history[mintime][statekey][item] - history[mintime + 1][statekey][item]) for item in inter])
           mysum += sum ([0.0] + [history[mintime][statekey][item] for item in dif1])
           mysum += sum ([0.0] + [history[mintime+1][statekey][item] for item in dif2])
       if mysum >= FRACEPSILON:
          return False 
    return True


def isStopCond(history,G,smodel,perfectinfo,typemethod):
    """should we stop inferring
       Args:
         history: all history inferred so far
         G: Graph
         smodel: spreading model
         perfectinfo: bool
         typemethod:
       Returns:
         flag: boolean: should algorithm stop or not
    """
    funcname = "is{0}StopCond".format(perfectinfo)
    return globals()[funcname](history,G,smodel,typemethod)


#INFTIME = 1000000
def frac2Trace(frachist,smodel):
    """converts fractional hist data to trace format
    Args:
       frachist:
       smodel:
    Returns:
       trace:
    """
    allnodes = set(node for time in frachist.keys() for state in frachist[time].keys() for node in frachist[time][state])
    alltimes = frachist.keys()
    trace = {node: {time : {} for time in alltimes} for node in allnodes}
    for time in alltimes:
        for state in frachist[time].keys():
            if smodel == "si" and state not in [Trace.SUSCEPTIBLE,Trace.INFECTED]:
               continue
            if smodel == "sir" and state not in [Trace.SUSCEPTIBLE,Trace.INFECTED,Trace.RECOVERED]:
               continue 
            for node in frachist[time][state]:
                trace[node][time][state] = frachist[time][state][node]
    return trace
