#Price Collecting Vertex Cover and Price Collecting Dominating Set Vertex Cover reconstruction code
import networkx as nx
import numpy as np
import scipy as sp
import scipy.optimize
import random
import math
import sys
sys.path.append("./lib")
import os
import time
from copy import deepcopy
import myutilities as myutil
from Trace import Trace
from Constraint import Constraint
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from HistoryEnsemble import DiscreteNoneEnsemble as DisEns
import HistoryUtil
import checkUtil
import DiscreteGreedy

class IterPCDSVCTest():
     #iter PCDSVC Test Class
     def __init__():
         return

     @staticmethod  
     def testConstraints(cons2item,item2cons,savecons2item,saveitem2cons,linear,quad,curstate,G,method): 
         """tests constraints
         Args:
            cons2item: modified
            item2cons: modified 
            savecons2item: real one
            saveitem2cons: real one
            linear,quad,curstate,G:
            method:
         Returns:
            bool: true/false
         """
         for itemset in savecons2item.values():
             for var in itemset:
                 assert linear.has_key(var) or quad.has_key(var)   
         if method == "cover":
            assert len(savecons2item.keys()) == len(linear.keys()) + len(quad.keys())
         elif method == "uncover":
             assert len(savecons2item.keys()) == len(linear.keys())    
         for cons in savecons2item.keys():
             if len(savecons2item[cons].intersection(set(quad.keys()))) == 0:
                seennodes = set()
                for item in savecons2item[cons]:
                    node = int(item.replace("x","").split("?")[0])
                    assert node in curstate[Trace.INFECTED] or node in curstate[Trace.RECOVERED]
                    seennodes.add(node)
                centernode = None    
                for node in seennodes:
                    if len(set(G.predecessors(node)+[node]).intersection(seennodes)) == len(seennodes):
                       centernode = node
                       break
                assert centernode != None    
             else:
                assert len(savecons2item[cons]) == 3
         assert len(cons2item.keys()) == 0 and len(item2cons.keys()) == 0 
         return True
     
     @staticmethod
     def testGreedySolution(cursol,curcons2item,curitem2cons,linear,quad,curstate,G,method):
         """tests greedy solution
         Args:
            cursol: 
            curcons2items: modified
            curitem2cons: modified
            linear:
            quad:
            curstate:
            G:
            method:
         Returns:
            bool: true/false  
         """
         curnodes = []
         dependG = nx.DiGraph(G.subgraph(list(curstate[Trace.INFECTED]) + list(curstate[Trace.RECOVERED])))
         for node in dependG.nodes():
             assert node in curnodes or len(set(dependG.predecessors(node)).intersection(curnodes)) > 0
           
     @staticmethod 
     def testGvcSolution(prestate,retvalues,G,curstate):
         """check GVC solution
         Args:
            prestate:
            retvalues:
            G:
            curstate:
         Returns:
            bool: true/false 
         """
         for val in retvalues.values():
             assert val in [-0.0,0.0,0.5,1.0]  
         singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
         doublevars = [var for var in retvalues.keys() if var.find("_") != -1]
         posinodes =  curstate[Trace.INFECTED].union(curstate[Trace.RECOVERED])
         zeronodes = posinodes.difference([int(singlevar.replace("x","").split("?")[0]) for singlevar in singlevars])
         onenodes =  set(int(singlevar.replace("x","").split("?")[0]) for singlevar in singlevars if retvalues[singlevar] == 1.0)
         halfnodes = set(int(singlevar.replace("x","").split("?")[0]) for singlevar in singlevars if retvalues[singlevar] == 0.5)
         assert len(zeronodes.union(onenodes).union(halfnodes)) == len(posinodes)
         assert sum([len(prestate[state]) for state in prestate.keys()]) == sum([len(curstate[state]) for state in curstate.keys()])
         assert len(set(curstate[Trace.SUSCEPTIBLE]).difference(set(prestate[Trace.SUSCEPTIBLE]))) == 0
         assert len(set(prestate[Trace.RECOVERED]).difference(set(curstate[Trace.RECOVERED]))) == 0
         return True
     
     @staticmethod
     def testCoefs(linear,quad,smodel,curstate):
         """test coefs
         Args:
            linear:
            quad:
            smodel:
            curstate:
         Returns:
            bool" true/false    
         """               
         for varname in linear.keys():
             state = varname.replace("x","").split("?")[1]
             assert state == Trace.INFECTED
             assert linear[varname] >= 0
             node = int(varname.replace("x","").split("?")[0])
             assert node not in curstate[Trace.SUSCEPTIBLE]
         if smodel == "sir":
            for var1,var2 in quad.keys(): 
                assert linear.has_key(var1) or linear.has_key(var2)
         allnodes = [int(varname.replace("x","").split("?")[0]) for varname in linear.keys()]
         assert len(allnodes) == len(set(allnodes))
         for var1,var2 in quad.keys():
             state1 = var1.replace("x","").split("?")[1]
             state2 = var2.replace("x","").split("?")[1]
             assert state1 == Trace.INFECTED and state1 == state2
             assert not quad.has_key((var2,var1))
             assert quad[(var1,var2)] >= 0
         return True
           

def genGvcMainConst(linear,quad,curstate,G):
    """generates GVC constraints
    Args:
       linear:
       quad:
       curstate:
       G:
    Returns:
       consstr: main constraint string
    """
    varcompare = lambda x1,x2: (x1,x2) if x1 <= x2 else (x2,x1)
    consstr = "Subject To \n"
    cindex = 0
    for var1,var2 in quad.keys():
        consstr += " cons{0}: 1.0 {1} + 1.0 {2} + 1.0 {1}_{2} >= 1.0 \n".format(cindex,var1,var2)
        cindex += 1 
    return consstr

            
def genGvcBoundConst(linear,quad):
    """generates GVC bound constraints
    Args:
       linear:
       quad:
    Returns:
       boundstr:
    """
    boundstr = "Bounds\n"
    boundstr += "\n".join(["0 <= {0} <= 1.0".format(varname) for varname in linear.keys()]) + "\n" 
    boundstr += "\n".join(["0 <= {0}_{1} <= 1.0".format(var1,var2) for var1,var2 in quad.keys()]) + "\n" 
    return boundstr
 

def genGvcObj(linear,quad):
    """generates GVC objective function 
    Args:
      linear:
      quad:
    Returns:
      objstr:
    """
    isPlus = lambda x: "+" if x >= 0 else " "
    objstr = "Minimize\n obj: "
    objstr += " ".join([" {0} {1} {2} ".format(isPlus(linear[var1]),linear[var1],var1) for var1 in linear.keys()])       
    objstr += " ".join([" {0} {1} {2}_{3} ".format(isPlus(quad[(var1,var2)]),quad[(var1,var2)],var1,var2) for var1,var2 in quad.keys()])
    return objstr


def getUnitProb(info):
    """
    """
    return 1.0 - math.exp(-1.0*info[1][0])

def estimateGvcCoefs(curstate,smodel,G):
    """estimate GVC coefs for program
    Args:
       curstate:
       smodel:
       G:
    Returns:
       linear:
       quad:
    """
    varcompare = lambda x1,x2: (x1,x2) if x1 <= x2 else (x2,x1)
    transdist = {"si":Trace.S2I,"sir":Trace.S2I}
    if smodel == "si":
       posistates = [Trace.INFECTED] 
    elif smodel in ["sir"]:
       posistates = [Trace.INFECTED,Trace.RECOVERED]
    usenodes = list(curstate[Trace.INFECTED]) + list(curstate[Trace.RECOVERED]) 
    linear,quad = {"x{0}?{1}".format(inode,Trace.INFECTED):0.0 for inode in usenodes},{} #quad (1-in) here
    for inode in curstate[Trace.INFECTED]: 
        varname = "x{0}?{1}".format(inode,Trace.INFECTED)  
        if smodel == "sir":
           linear.setdefault(varname,0.0)
           linear[varname] += -1.0 * math.log(1.0 - getUnitProb(G.node[inode][Trace.I2R]))
        for neighnode in G.predecessors(inode):
            for state in posistates:
                if neighnode in curstate[state]:
                   varname2 = "x{0}?{1}".format(neighnode,Trace.INFECTED)
                   quadkey = varcompare(varname,varname2)      
                   quad.setdefault(quadkey,0.0)
                   quad[quadkey] += -1.0 * math.log(1.0 - getUnitProb(G[neighnode][inode][transdist[smodel]]))
    for rnode in curstate[Trace.RECOVERED]:
        varname = "x{0}?{1}".format(rnode,Trace.INFECTED)
        linear.setdefault(varname,0.0)
        linear[varname] += -1.0 * math.log(getUnitProb(G.node[rnode][Trace.I2R]))
    for snode in curstate[Trace.SUSCEPTIBLE]:
        for neighnode in G.predecessors(snode):
            for state in posistates:
                if neighnode in curstate[state]:
                   varname = "x{0}?{1}".format(neighnode,Trace.INFECTED)
                   linear.setdefault(varname,0.0)
                   linear[varname] += -1.0 * math.log(1.0 - getUnitProb(G[neighnode][snode][transdist[smodel]]))     
    return linear,quad              


def runGvc(curstate,smodel,G,runfolder,index,restrict): 
    """runs GVC at each time step
    Args:
       curstate:
       smodel:
       G:
       runfolder:
       index:
       restrict:
    Returns:
       sol:
    """ 
    linear,quad = estimateGvcCoefs(curstate,smodel,G)
    objstr = genGvcObj(linear,quad)
    consstr = genGvcMainConst(linear,quad,curstate,G)
    boundstr = genGvcBoundConst(linear,quad)
    print "prerun info", len(curstate[Trace.SUSCEPTIBLE]),len(curstate[Trace.INFECTED]),len(curstate[Trace.RECOVERED])
    retvalues,objval,runtime = HistoryUtil.runHistoryCode(consstr,objstr,boundstr,runfolder,["x"],"pcvc-{0}".format(index),"removeZeros")
    sol = roundGvcSolution(retvalues,curstate,G)
    if globals()["isTest"]:
       IterPCDSVCTest.testCoefs(linear,quad,smodel,curstate) 
       IterPCDSVCTest.testGvcSolution(sol,retvalues,G,curstate) 
       checkUtil.isSnapshotValid(sol,curstate,G,smodel,"notdetailed") 
    return sol


def roundGvcSolution(retvalues,curstate,G):
    """round randomized solution to gvc solution
    Args:
       retvalues:
       curstate:
       G:
    Returns:
       prestate:
    """  
    prestate = {Trace.SUSCEPTIBLE: deepcopy(curstate[Trace.SUSCEPTIBLE]),Trace.EXPOSED:set(),Trace.INFECTED:set(),Trace.RECOVERED:set()}        
    singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
    posinodes =  curstate[Trace.INFECTED].union(curstate[Trace.RECOVERED])
    zeronodes = posinodes.difference([int(singlevar.replace("x","").split("?")[0]) for singlevar in singlevars])
    onenodes =  set(int(singlevar.replace("x","").split("?")[0]) for singlevar in singlevars if retvalues[singlevar] == 1.0)
    halfnodes = set(int(singlevar.replace("x","").split("?")[0]) for singlevar in singlevars if retvalues[singlevar] == 0.5)
    for halfnode in halfnodes:
        if random.random() < 0.5:
           onenodes.add(halfnode)
        else:
           zeronodes.add(halfnode)     
    prestate[Trace.INFECTED] |= set(onenodes)
    for zeronode in zeronodes:
        if zeronode in curstate[Trace.INFECTED]:
           prestate[Trace.SUSCEPTIBLE].add(zeronode)
        elif zeronode in curstate[Trace.RECOVERED]:
           prestate[Trace.RECOVERED].add(zeronode)
    return prestate


def prepareConstraint(linear,quad,curstate,G,algo):
    """prepares constraints
    Args: 
       linear:
       quad: 
       curstate:
       G:
       algo:
    Returns:
       cons2item:
       item2cons:
    """    
    cons2item = {}
    cons = 0
    for var1,var2 in quad.keys():
        cons2item[cons] = set([(var1,var2),var1,var2])
        cons += 1
    if algo == "pcdsvc":
       posnodes = curstate[Trace.INFECTED].union(curstate[Trace.RECOVERED])
       for node in posnodes:
           cons2item[cons] = set("x{0}?{1}".format(subnode,Trace.INFECTED) for subnode in set(G.predecessors(node)+[node]).intersection(posnodes))
           cons += 1
    item2cons = {}      
    for cons in cons2item.keys():
        for item in cons2item[cons]:
            item2cons.setdefault(item,set())
            item2cons[item].add(cons)
    return cons2item,item2cons

               
def runGreedy(curstate,smodel,G,runfolder,algo,restrict):
    """greedy PCDSVC or pcvc reconstruction
    Args:
       curstate:   
       smodel:
       G:
       runfolder:
       algo:
       restrict:
    Returns:
       sol:
    """    
    linear,quad = estimateGvcCoefs(curstate,smodel,G)
    cons2item, item2cons = prepareConstraint(linear,quad,curstate,G,algo)
    if globals()["isTest"]:
       savecons2item, saveitem2cons = deepcopy(cons2item), deepcopy(item2cons)     
    print "mymyinfo"    
    print len(cons2item.keys())
    print len(item2cons.keys())
    print len(linear.keys())
    print len(quad.keys())  
    varmap = {}
    for var in linear.keys():
        node,state = var.replace("x","").split("?")
        #print var,node,state
        node = int(node)
        assert not varmap.has_key(node)
        varmap[node] = state
    print len(varmap.keys())
    print len(set(list(curstate[Trace.INFECTED])+list(curstate[Trace.RECOVERED])))    
    assert len(set(varmap.keys()).symmetric_difference(set(list(curstate[Trace.INFECTED])+list(curstate[Trace.RECOVERED])))) == 0    
    remforbidX = set() #items that can not be removed
    addforbidX = set() #items that can not be added 
    for node in set(restrict.keys()).intersection(set(varmap.keys())):
        if varmap[node] == restrict[node]:
           varname = "x{0}?{1}".format(node,restrict[node])  
           remforbidX.add(varname)
           if item2cons.has_key(varname):
              mycons = item2cons[varname]
              for cons in mycons:
                  for item in cons2item[cons]:
                      item2cons[item].remove(cons)
              for cons in mycons:
                  del cons2item[cons]              
              del item2cons[varname]
           for item in item2cons.keys():
               if len(item2cons[item]) == 0:
                  del item2cons[item] 
        else:
           varname = "x{0}?{1}".format(node,varmap[node]) 
           addforbidX.add(varname)     
    print remforbidX
    print addforbidX
    print len(remforbidX)
    print len(addforbidX)
    curobj = 0.0   
    cursol = []
    if len(remforbidX) != 0:
       cursol = list(remfordbidX)
       curobj = sum([linear[varname] for varname in cursol if linear.has_key(varname)]) 
    while True:
        minitem = None
        mincoef = 100000000.0 
        minitemobj = None
        for item in set(item2cons.keys()).difference(addforbidX):
            if len(item2cons[item]) == 0:
               continue 
            if linear.has_key(item):
               difobj = linear[item]
            elif quad.has_key(item):
               difobj = quad[item] 
            itemcoef = difobj / float(len(item2cons[item]))
            if itemcoef < mincoef:
               mincoef = itemcoef
               minitem = item
               minitemobj = curobj + difobj
        assert minitem != None
        curobj = minitemobj
        newcons = list(item2cons[minitem])
        for cons in newcons:
            for titem in cons2item[cons]:
                item2cons[titem].remove(cons)
        del item2cons[minitem]
        for cons in newcons:
            del cons2item[cons]
        for item in item2cons.keys():
            if len(item2cons[item]) == 0:
               del item2cons[item]    
        cursol.append(minitem)  
        if len(cons2item.keys()) == 0:
           break
    print len(cursol)
    print cursol
    edges = []
    nodes = []
    if globals()["isTest"]:
       IterPCDSVCTest.testConstraints(cons2item,item2cons,savecons2item,saveitem2cons,linear,quad,curstate,G,algo) 
       #IterPCDSVCTest.testGreedySolution(cursol,cons2item,item2cons,linear,quad,curstate,method)
       IterPCDSVCTest.testCoefs(linear,quad,smodel,curstate) 
       checkUtil.isSnapshotValid(cursol,curstate,G,smodel,"notdetailed") 
    print cursol  
    return cursol,curobj


def inferPreSnapshot(algo,G,smodel,algoinfo,seenstate,runfolder,infermode):
    """infers pre snapshot before t_{min}
    Args:
       algo: 
       G:      
       smodel:
       algoinfo:
       seenstate:
       runfolder:
       infermode:
    Returns:
       history:
    """
    enscount = 1
    ensemble = algoinfo["ensemble"]
    if ensemble:
       enscount = algoinfo["enscount"]   
    allsol = []
    for index in xrange(enscount):
        curstate = deepcopy(seenstate)
        sol,curcount= {}, 0
        while not HistoryUtil.isStopCond(sol,G,smodel,"Perfect",algo):
           if algo == "Pcvc" and globals()["ALGO"] == "lp":
              curstate,tobj = runGvc(curstate,smodel,G,runfolder,curcount,{})
           elif algo == "Pcvc" and globals()["ALGO"] == "greedy":
              curstate,tobj = runGreedy(curstate,smodel,G,runfolder,algo,{})  
           elif algo == "Pcdsvc":
              curstate,tobj = runGreedy(curstate,smodel,G,runfolder,algo,{})       
           sol[-1 * (curcount + 1)] = deepcopy(curstate)
           curcount += 1
        allsusflag = True    
        for time in sorted(sol.keys(),reverse=True):
            if len(sol[time][Trace.SUSCEPTIBLE]) != G.number_of_nodes():
               allsusflag = False
               break
        if allsusflag:
           sol = {-1:deepcopy(seenstate)}
           allsol.append(sol)
           continue 
        print "sol"   
        for time in sorted(sol.keys()):
            print time,len(sol[time][Trace.SUSCEPTIBLE]),len(sol[time][Trace.EXPOSED]),len(sol[time][Trace.INFECTED]),len(sol[time][Trace.RECOVERED])
        sol = HistoryUtil.lastFixSolution(sol,G,smodel,"None",infermode)
        print "after"
        for time in sorted(sol.keys()):
            print time,len(sol[time][Trace.SUSCEPTIBLE]),len(sol[time][Trace.EXPOSED]),len(sol[time][Trace.INFECTED]),len(sol[time][Trace.RECOVERED])
        for time in sol.keys():
            for state in allstates:
                sol[time].setdefault(state,set())
        allsol.append(sol)
    if globals()["isTest"]:    
       for sol in allsol:
           checkUtil.checkPreConversion(sol,G,smodel) 
    if ensemble:
       return allsol 
    elif infermode == "Spreader":
       return allsol[0][min(allsol[0].keys())], min(allsol[0].keys())
    elif infermode == "History":
       return allsol[0]

    
def inferBiDirection(prestate,nextstate,smodel,G,restrict):
    """infers bidirectional
    Args:
       prestate:
       nextstate:
       smodel:
       G:
       restrict:
    Returns:
    """
    return
    
    
def inferForward(curstate,smodel,G,restrict):
    """infers forward way
    Args:
       curstate:
       smodel:
       G:
       restrict:
    Returns:
       inferstate: state at the next step
       inferscore:
    """
    return DiscreteGreedy.inferForward(curstate,smodel,G,restrict)


def getRestrictions(lowtstate,hightstate,smodel):
    """get restrictions given higher and lower bounds
    Args:
       hightstate: high t
       lowtstate: low t
       smodel:
    Returns:
       restrict:   
    """
    return DiscreteGreedy.getRestrictions(lowtstate,hightstate,smodel)
 

def inferBetweenDP(prestate,nextstate,pretime,nexttime,G,smodel,algo,algoinfo):
    """infer the history between consecutive snapshots by dynamic programming
    Args:
       prestate:
       nextstate:
       pretime:
       nexttime:
       G:
       smodel:
       algo:
       algoinfo:
    Returns:
       inferhist:
    """
    assert nexttime-pretime >= 2
    print "info",pretime,nexttime
    forhist = {time:None for time in xrange(pretime,nexttime-1)}
    forhist[pretime] = deepcopy(prestate)
    backhist = {time:None for time in xrange(pretime+2,nexttime+1)}
    backhist[nexttime] = deepcopy(nextstate)
    forhistscore = {time:None for time in xrange(pretime,nexttime-1)}
    backhistscore = {time:None for time in xrange(pretime+2,nexttime+1)}
    forhistscore[pretime] = 0.0
    backhistscore[nexttime] = 0.0
    print "forward:"
    for fortime in xrange(pretime+1,nexttime-1):  
        assert forhist[fortime-1] != None and forhist[fortime] == None
        restrict = getRestrictions(forhist[fortime-1],nextstate,smodel)
        forhist[fortime],forscore = [deepcopy(item) for item in inferForward(forhist[fortime-1],smodel,G,restrict)]
        forhistscore[fortime] = -1.0*forscore
    print "backward:"    
    for backtime in xrange(nexttime-1,pretime+1,-1):    
        assert backhist[backtime] == None and backhist[backtime+1] != None
        objclass = GreedyObjFunc(G,smodel,"backward",[backhist[backtime+1]])
        objmethod = getattr(objclass,"get{0}ObjFunc".format(algoinfo["objmethod"].capitalize()))
        restrict = getRestrictions(prestate,backhist[backtime+1],smodel)
        consts = genCoverCons(G,smodel,backhist[backtime+1])
        if algo == "Pcvc" and globals()["ALGO"] == "lp":
           tstate,backscore = runGvc(curstate,smodel,G,runfolder,curcount,restrict)
        elif algo == "Pcvc" and globals()["ALGO"] == "greedy":
           tstate,backscore = runGreedy(curstate,smodel,G,runfolder,algo,restrict)  
        elif algo == "Pcdsvc":
           tstate,backscore = runGreedy(curstate,smodel,G,runfolder,algo,restrict) 
        backhist[backtime] = deepcopy(tstate)
        backhistscore[backtime] = backscore
    maxscore,maxbipoint,maxbistate = None,None,None
    for bipoint in xrange(pretime+1,nexttime):
        assert forhist[bipoint-1] != None and backhist[bipoint+1] != None and forhistscore[bipoint-1] != None and backhistscore[bipoint+1] != None
        print "bifind",bipoint,pretime,nexttime
        if len(backhist[bipoint+1][Trace.SUSCEPTIBLE].difference(forhist[bipoint-1][Trace.SUSCEPTIBLE])) != 0:
           continue
        restrict = getRestrictions(forhist[bipoint-1],backhist[bipoint+1],smodel)
        tempbistate,biscore = [deepcopy(item) for item in inferBiDirection(forhist[fortime-1],smodel,G,restrict)]
        curscore = forhistscore[bipoint-1] + backhistscore[bipoint+1] + biscore
        if maxscore == None or curscore >= maxscore:
           maxscore = curscore 
           maxbipoint = bipoint
           maxbistate = deepcopy(tempbistate)
    assert maxbipoint != None and maxbistate != None   
    inferhist = {maxbipoint: deepcopy(maxbistate)}
    for curtime in xrange(pretime+1,maxbipoint):
        assert not inferhist.has_key(curtime)
        inferhist[curtime] = deepcopy(forhist[curtime])
    for curtime in xrange(maxbipoint+1,nexttime):
        assert not inferhist.has_key(curtime)
        inferhist[curtime] = deepcopy(backhist[curtime])
    if globals()["isTest"]:    
       assert DiscreteGreedyTest.testForwardInference(forhist,prestate,nextstate)    
       assert len(set(range(pretime+1,nexttime)).symmetric_difference(set(inferhist.keys()))) == 0    
    return inferhist


isTest = True
ALGO = "lp" #"greedy"
def runner(algo,algoinfo,curstates,G,smodel,infermode,resultfile,runfolder,obstimes,TESTMODE=False):
    """iterative pcdsvc and pcvc relaxations
    Args:
       algo:
       algoinfo:
       curstates: current state knowledge
       G: Graph
       smodel: spreading model
       infermode: infermode
       resultfile:
       runfolder:
       TESTMODE: 
    Returns:
       sol: solution
    """
    assert smodel in ["si","sir"] and checkUtil.checkPreTrace(curstates,obstimes) and sorted(obstimes) == obstimes and len(obstimes) == len(curstates) and obstimes[0] == 0
    globals()["isTest"] = TESTMODE
    assert checkUtil.checkPreTrace(curstates,obstimes)
    temphistory = inferPreSnapshot(algo,G,smodel,algoinfo,deepcopy(curstates[0]),runfolder,infermode)
    if infermode == "Spreader":
       if algoinfo["ensemble"]: 
          enscount = algoinfo["enscount"]
          assert enscount == len(temphistory)     
          history = deepcopy(DisEns.getEnsembleSolution(temphistory,infermode,curstates,obstimes,G,smodel))
       else:   
          history = deepcopy(temphistory)
    elif infermode == "History": 
       for state in curstates[0].keys():
           assert len(curstates[0][state]) == len(temphistory[0][state]) 
       assert not temphistory.has_key(0)     
       history = {}
       for time in temphistory.keys():
           history[time+inittime] = deepcopy(temphistory[time]) 
       for index in xrange(1,len(obstimes)):
           pretime,nexttime = obstimes[index-1:index+1]
           prestate,nextstate = curstates[index-1:index+1]
           if nexttime-pretime == 1: 
              continue 
           print "iter: ",index,pretime,nexttime
           temphistory = inferBetweenDP(deepcopy(prestate),deepcopy(nextstate),pretime,nexttime,G,smodel,algo,algoinfo)
           for time in temphistory.keys():
               assert not history.has_key(time+inittime)
               history[time+inittime] = deepcopy(temphistory[time])
       for tindex in xrange(len(obstimes)):
           time = obstimes[tindex]
           assert not history.has_key(time+inittime)
           history[time+inittime] = deepcopy(curstates[tindex]) 
       assert len(set(history.keys()).symmetric_difference(range(0,max(history.keys())+1))) == 0    
       print len(history.keys())
       print "print hist"              
       for time in sorted(history.keys()):
           for state in history[time].keys():
               print time,state,len(history[time][state])
       print inittime
       print obstimes        
       history = InputOutput.history2Times(history,curstates,smodel)
    if globals()["isTest"]:
       checkUtil.checkSolution(infermode,"None",history,G,smodel,"notdetailed")  
    InputOutput.writeHistoryResult(history,infermode,resultfile)    
    
 
if __name__ == "__main__":
    with gzip.open(sys.argv[1],"rb") as infile: 
        algo = cPickle.load(infile)
        algoinfo = cPickle.load(infile)       
        curstates = cPickle.load(infile)
        G = cPickle.load(infile)
        smodel = cPickle.load(infile)
        infermode = cPickle.load(infile)
        resultfile = cPickle.load(infile)
        runfolder = cPickle.load(infile)
        obstimes = cPickle.load(infile)
    runner(algoinfo,curstates,G,smodel,infermode,resultfile,runfolder)
