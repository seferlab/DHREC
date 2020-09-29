#Iterative Min Cut history relaxation reconstruction code
#Suitable for seir model
#probability of 0 infected nodes(we must modify while loop accordingly) ?
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
import numpy as np
import checkUtil
import cutUtil



class IterMinCutTest():
     #iter min cut Test Class
     def __init__():
         return

     @staticmethod          
     def checkLpMinCutSolution(prestate,retvalues,G,smodel,curstate):
         """check mincut solution
         Args:
            prestate:
            retvalues:
            G:
            smodel:
            curstate:
         Returns:
            bool: true/false 
         """
         checkUtil.isSnapshotValid(prestate,curstate,G,smodel,"notdetailed")    
         singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
         doublevars = [var for var in retvalues.keys() if var.find("_") != -1]
         seennodes = []
         for singlevar in singlevars:
             assert retvalues[singlevar] == 1.0
             node = int(singlevar.replace("x","").split("?")[0])
             seennodes.append(node)
         assert len(set(seennodes)) == len(seennodes) 
         return True

     @staticmethod
     def checkFlowMinCutSolution(prestate,G,curstate,smodel):
         """check mincut solution
         Args:
           prestate:
           G:
           curstate:
           smodel:
         Returns:
           bool: true/false 
         """
         checkUtil.isSnapshotValid(prestate,curstate,G,smodel,"notdetailed")
         return True

     @staticmethod
     def AllSolutionCheck(flowsol,nxflowsol,lpsol,flow,nxflow,curstate,smodel,cutG):
         """checks all solutions
         Args:
            flowsol:
            nxflowsol:
            lpsol:  
            flow:
            nxflow:
            curstate:
            smodel:
            cutG:
         Returns:
            bool: true/false
         """
         flowtransnodes = HistoryUtil.assignStateTransitions(curstate,flowsol,smodel)
         nxflowtransnodes = HistoryUtil.assignStateTransitions(curstate,nxflowsol,smodel)
         lptransnodes = HistoryUtil.assignStateTransitions(curstate,lpsol,smodel)
        
         s,t = min(cutG.nodes()), max(cutG.nodes())
         snodes,tnodes = set([s]), set([t])
         lpcut = 0.0  
         for node in lptransnodes[(Trace.SUSCEPTIBLE,Trace.EXPOSED)]:
             tnodes.add(node)
         for node in lptransnodes[(Trace.EXPOSED,Trace.EXPOSED)]:
             snodes.add(node)         
         for node in lptransnodes[(Trace.EXPOSED,Trace.INFECTED)]:
             snodes.add(node)
         for node in lptransnodes[(Trace.INFECTED,Trace.INFECTED)]:
             tnodes.add(node)
         for node in lptransnodes[(Trace.INFECTED,Trace.RECOVERED)]:
             tnodes.add(node)
         for node in lptransnodes[(Trace.RECOVERED,Trace.RECOVERED)]:
             snodes.add(node)
         assert len(snodes.intersection(tnodes)) == 0 and len(snodes) + len(tnodes) == cutG.number_of_nodes()
         for snode in snodes:
             for tnode in tnodes:
                 if cutG.has_edge(snode,tnode):
                    lpcut += cutG[snode][tnode]["weight"]
         assert abs(lpcut-flow) < 0.1 and abs(lpcut-nxflow) < 0.1 and abs(nxflow-flow) < 0.1        
         return True

     
def genMinCutMainConst(linear,quad,smodel,G):
    """generates MinCut constraints
    Args:
       linear:
       quad:
       smodel:
       G:
    Returns:
       consstr: main constraint string
    """
    consstr = "Subject To \n"
    for var1,var2 in quad.keys():
        consstr += " {0} - {1} - {1}_{0} <= 0.0 \n".format(var2,var1)
    return consstr
     
    inodes = []
    if smodel == "seir":
       for varname in linear.keys(): 
           splitted = varname.replace("x","").split("?")
           if splitted[1] == Trace.INFECTED:
              inodes.append(int(splitted[0]))
    for varname in linear.keys():
        [node,state] = varname.replace("x","").split("?")
        if state != Trace.SUSCEPTIBLE:
           continue
        node = int(node)
        linestr = " 1.0 {0} <= ".format(varname)
        strlist = []
        for neighnode in set(G.predecessors(node)).intersection(set(inodes)):
            neighvar = "x{0}?{1}".format(neighnode,Trace.INFECTED)
            strlist.append(" 1.0 {0} ".format(neighvar))
        if len(strlist) != 0:     
           linestr += " + ".join(strlist) + " \n "
        else:
           linestr += " 0.0 \n "
        consstr += linestr 
    return consstr

        
def genMinCutBoundConst(linear,quad):
    """generates MinCut bound constraints
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


def genMinCutObj(linear,quad):
    """generates MinCut objective function 
    Args:
      linear:
      quad:
    Returns:
      objstr:
    """
    isPlus = lambda x: "+" if x >= 0 else " "
    objstr = "Minimize\n obj: "
    objstr += " ".join([" {0} {1} {2} ".format(isPlus(linear[var1]),linear[var1],var1) for var1 in linear.keys()])       
    objstr += " ".join([" {0} {1} {2} ".format(isPlus(quad[(var1,var2)]),quad[(var1,var2)],var1+"_"+var2) for var1,var2 in quad.keys()])
    return objstr


def getUnitProb(info):
    """
    """
    return 1.0 - math.exp(-1.0*info[1][0])


def getMinCutFromFlow(cutG,F,s,t):
    """returns mincut from flow graph
    Args:
      cutG: graph
      F: flow info
      s: source 
      t: destination
    Returns:
      node2label:
    """
    temp = nx.DiGraph(cutG)
    for node1,node2 in cutG.edges():
        if cutG[node1][node2]["weight"] == F[node1][node2]:
           temp.remove_edge(node1,node2)
    node2label = {node:t for node in cutG.nodes()}        
    curset,seenset = [s],set([s])
    while True:
       if len(curset) == 0:
          break
       curitem = curset.pop() 
       difset = set(temp.successors(curitem)).difference(seenset)
       seenset |= difset
       curset.extend(difset)
    for node in seenset:
        node2label[node] = s
        
    if globals()["isTest"]:
       assert cutG.number_of_nodes() == len(node2label.keys())
    return node2label


def runBoykovFlow(cutG,G,curstate,runfolder):
    """runs networkx flow algo
    Args:
       cutG:
       G:
       curstate:
       runfolder:
    Returns:
       sol:
    """
    s,t = min(cutG.nodes()), max(cutG.nodes()) 
    CODEPATH = "{0}/{1}".format(globals()["BOYKOVCODEDIR"],globals()["BOYKOVCODENAME"])
    mymap ={"s":s,"t":t}
    inputfile = "{0}/data{1}.input".format(runfolder,random.random())
    outfile = "{0}/data{1}.out".format(runfolder,random.random())    
    cutUtil.storeBoykovGraph(cutG,s,t,inputfile)
    code = "{0} {1} {2}".format(CODEPATH,inputfile,outfile)
    os.system(code)
    node2label,count,flow = {},0, None
    with open(outfile,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if count == 0:
               flow = float(line) 
               count +=1
               continue
            node,terminal = line.split("\t")
            node2label[int(node)] = mymap[terminal] 
    for filename in [inputfile,outfile]:
        os.system("rm -rf {0}".format(filename)) 
    prestate = {Trace.SUSCEPTIBLE:deepcopy(curstate[Trace.SUSCEPTIBLE]),Trace.EXPOSED:set(),Trace.INFECTED:set(),Trace.RECOVERED:set()}
    mapback = {Trace.EXPOSED: {s:Trace.EXPOSED, t:Trace.SUSCEPTIBLE},Trace.INFECTED:{s:Trace.EXPOSED, t:Trace.INFECTED},Trace.RECOVERED:{s:Trace.RECOVERED, t:Trace.INFECTED}}
    for node in set(node2label.keys()).difference(set([s,t])):
        for state in curstate.keys():
            if node in curstate[state]:
               prestate[mapback[state][node2label[node]]].add(node)
               break
    return flow,prestate 


def runNetworkxFlow(cutG,G,curstate):
    """runs flow algo
    Args:
       cutG:
       G:
       curstate:
    Returns:
       sol:
    """
    s,t = min(cutG.nodes()), max(cutG.nodes())
    flow, F = nx.ford_fulkerson(cutG,s,t,capacity="weight")
    node2label = getMinCutFromFlow(cutG,F,s,t)
    prestate = {Trace.SUSCEPTIBLE:deepcopy(curstate[Trace.SUSCEPTIBLE]),Trace.EXPOSED:set(),Trace.INFECTED:set(),Trace.RECOVERED:set()}
    mapback = {Trace.EXPOSED: {s:Trace.EXPOSED, t:Trace.SUSCEPTIBLE},Trace.INFECTED:{s:Trace.EXPOSED, t:Trace.INFECTED},Trace.RECOVERED:{s:Trace.RECOVERED, t:Trace.INFECTED}}
    for node in set(node2label.keys()).difference(set([s,t])):
        for state in curstate.keys():
            if node in curstate[state]:
               prestate[mapback[state][node2label[node]]].add(node)
               break
    return flow,prestate 
    
    
def buildCutGraph(linear,quad,G):
    """builds min cut graph
    Args:
       linear:
       quad:
       G:
    Returns:
       cutG: cut graph 
    """
    cutG = nx.DiGraph()
    for var in linear.keys():
        node = int(var.replace("x","").split("?")[0])
        cutG.add_node(node)
    snode, tnode = min(G.nodes()) -1, max(G.nodes()) +1
    cutG.add_node(snode)
    cutG.add_node(tnode)   
    for var in linear.keys():
        node = int(var.replace("x","").split("?")[0])
        if linear[var] >= 0.0:
           cutG.add_edge(snode,node,weight=linear[var])
        else:
           cutG.add_edge(node,tnode,weight=-1.0*linear[var]) 
    for var1,var2 in quad.keys():
        node1 = int(var1.replace("x","").split("?")[0])
        node2 = int(var2.replace("x","").split("?")[0])
        assert not cutG.has_edge(node1,node2) and not cutG.has_edge(node2,node1) 
        cutG.add_edge(node2,node1,weight=quad[(var1,var2)])
        if cutG.has_edge(snode,node2):  
           cutG[snode][node2]["weight"] += quad[(var1,var2)]  
        else:    
           cutG.add_edge(snode,node2,weight=quad[(var1,var2)]) 
        if cutG.has_edge(node1,tnode): 
           cutG[node1][tnode]["weight"] += quad[(var1,var2)] 
        else:   
           cutG.add_edge(node1,tnode,weight=quad[(var1,var2)])  
    if globals()["isTest"]:
       curvars = {Trace.SUSCEPTIBLE:set(), Trace.INFECTED:set()}
       for varname in linear.keys():
           node, state = varname.replace("x","").split("?")
           node = int(node)
           curvars[state].add(node)
       for sendnode,recnode in cutG.edges():
           if sendnode == snode or recnode == tnode:
              continue
           assert sendnode in curvars[Trace.SUSCEPTIBLE] and recnode in curvars[Trace.INFECTED]      
    return cutG


def estimateCoefs(curstate,smodel,G):
    """estimate coefs for program
    Args:
       curstate:
       smodel:
       G:
    Returns:
       linear:
       quad:
    """
    assert smodel in ["seir"]
    linear,quad = {},{} #quad s(1-i)
    statemap = {Trace.INFECTED:Trace.INFECTED,Trace.RECOVERED:Trace.INFECTED}
    for enode in curstate[Trace.EXPOSED]: 
        varname = "x{0}?{1}".format(enode,Trace.SUSCEPTIBLE)
        linear.setdefault(varname,0.0)
        linear[varname] += math.log(1.0 - getUnitProb(G.node[enode][Trace.E2I]))
        for neighnode in G.predecessors(enode):
            for state in statemap.keys():
                if neighnode in curstate[state]:
                   varname2 = "x{0}?{1}".format(neighnode,statemap[state])
                   quadkey = (varname2,varname)
                   quad.setdefault(quadkey,0.0)
                   quad[quadkey] += -1.0 * math.log(1.0 - getUnitProb(G[neighnode][enode][Trace.S2E]))              
    for inode in curstate[Trace.INFECTED]: 
        varname = "x{0}?{1}".format(inode,Trace.INFECTED)
        linear.setdefault(varname,0.0)
        mydiv = getUnitProb(G.node[inode][Trace.E2I]) / (1.0 - getUnitProb(G.node[inode][Trace.I2R]))
        linear[varname] += math.log(mydiv)
    for rnode in curstate[Trace.RECOVERED]:
        varname = "x{0}?{1}".format(rnode,Trace.INFECTED)
        linear.setdefault(varname,0.0)
        linear[varname] += -1.0 * math.log(getUnitProb(G.node[rnode][Trace.I2R]))
    for snode in curstate[Trace.SUSCEPTIBLE]:
        for neighnode in G.predecessors(snode):
            for state in statemap.keys():
                if neighnode in curstate[state]:
                   varname = "x{0}?{1}".format(neighnode,statemap[state])
                   linear.setdefault(varname,0.0)
                   linear[varname] += -1.0 * math.log(1.0 - getUnitProb(G[neighnode][snode][Trace.S2E]))     
    if globals()["isTest"]:
       allnodes = [int(varname.replace("x","").split("?")[0]) for varname in linear.keys()]
       assert len(allnodes) == len(set(allnodes))
       for var1,var2 in quad.keys():
           assert quad[(var1,var2)] >= 0.0
           assert not quad.has_key((var2,var1))
           node1 = int(var1.replace("x","").split("?")[0])
           node2 = int(var2.replace("x","").split("?")[0])
           assert node2 in curstate[Trace.EXPOSED] and (node1 in curstate[Trace.INFECTED] or node1 in curstate[Trace.RECOVERED])
           assert linear.has_key(var1) and linear.has_key(var2) 
    return linear,quad


  
RUNMODE = "test"  #"flow","lp","flownetworkx", "test"
BOYKOVCODEDIR = "maxflow-v3.02"
BOYKOVCODENAME = "boykovcut"
def runMinCut(curstate,smodel,G,runfolder): 
    """runs min cut at each time step
    Args:
       curstate:
       smodel:
       G:
       runfolder:
    Returns:
       sol:
    """
    assert RUNMODE in ["flow","lp","flownetworkx","test"]
    if len(curstate[Trace.SUSCEPTIBLE]) == G.number_of_nodes():
       return curstate
    linear,quad = estimateCoefs(curstate,smodel,G)
    if globals()["RUNMODE"] in ["flow","test"]:
       cutG = buildCutGraph(linear,quad,G)
       flow,flowsol = runBoykovFlow(cutG,G,curstate,runfolder)
       if globals()["isTest"]:
          IterMinCutTest.checkFlowMinCutSolution(flowsol,G,curstate,smodel) 
    if globals()["RUNMODE"] in ["flownetworkx","test"]:
       cutG = buildCutGraph(linear,quad,G)
       nxflow,nxflowsol = runNetworkxFlow(cutG,G,curstate)
       if globals()["isTest"]:
          IterMinCutTest.checkFlowMinCutSolution(nxflowsol,G,curstate,smodel) 
    if globals()["RUNMODE"] in ["lp","test"]:
       objstr = genMinCutObj(linear,quad)
       consstr = genMinCutMainConst(linear,quad,smodel,G)
       boundstr = genMinCutBoundConst(linear,quad)
       retvalues,objval,runtime = HistoryUtil.runHistoryCode(consstr,objstr,boundstr,runfolder,["x"],"mincut","removeZeros")
       lpsol = convertLpSolution(retvalues,curstate)
       if globals()["isTest"]:
          IterMinCutTest.checkLpMinCutSolution(lpsol,retvalues,G,smodel,curstate)
          
    if globals()["RUNMODE"] == "test":
       IterMinCutTest.AllSolutionCheck(flowsol,nxflowsol,lpsol,flow,nxflow,curstate,smodel,cutG)
    if globals()["RUNMODE"] == "flownetworkx":
       return nxflowsol
    elif globals()["RUNMODE"] == "lp":
       return lpsol 
    elif globals()["RUNMODE"] in ["flow","test"]:
       return flowsol 


def convertLpSolution(retvalues,curstate):
    """converts solution
    Args:
       retvalues:
       curstate:
    Returns:
       prestate:
    """          
    singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
    if globals()["isTest"]:
       doublevars = [var for var in retvalues.keys() if var.find("_") != -1]
       for singlevar in singlevars:
           assert retvalues[singlevar] == 1.0 
       for doublevar in doublevars:
           assert retvalues[doublevar] == 1.0 
    prestate = {Trace.SUSCEPTIBLE: deepcopy(curstate[Trace.SUSCEPTIBLE]),Trace.EXPOSED:set(),Trace.INFECTED:set(),Trace.RECOVERED:set()} 
    for singlevar in singlevars:
        splitted = singlevar.replace("x","").split("?")
        node,state = int(splitted[0]),splitted[1]
        prestate[state].add(node)
    remainnodes = [node for state in [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED] for node in curstate[state]]                
    for node in remainnodes:
        if node in curstate[Trace.EXPOSED] and node not in prestate[Trace.SUSCEPTIBLE]:
           prestate[Trace.EXPOSED].add(node)
        elif node in curstate[Trace.INFECTED] and node not in prestate[Trace.INFECTED]:
           prestate[Trace.EXPOSED].add(node) 
        elif node in curstate[Trace.RECOVERED] and node not in prestate[Trace.INFECTED]:
           prestate[Trace.RECOVERED].add(node) 
    if globals()["isTest"]:
       tnodes = set(node for state in prestate.keys() for node in prestate[state])
       tnodes2 = set(node for state in curstate.keys() for node in curstate[state])
       assert len(tnodes) == len(tnodes2)
    return prestate


def inferBidirection(curstate,smodel,G,restrict):
    """infers both way
    Args:
       
    Returns:
    """
    return


def inferForward(curstate,smodel,G,restrict):
    """infers in forward direction
    Args:
       curstate: 
       smodel: 
       G:
       restrict:
    Returns:
       inferstate: state at the next step
       inferscore:
    """
    getCoefs = {}
    return


def inferBetweenDP(prestate,nextstate,pretime,nexttime,G,smodel,algoinfo):
    """infer the history between consecutive snapshots by dynamic programming
    Args:
       prestate:
       nextstate:
       pretime:
       nexttime:
       G:
       smodel:
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
        forhistscore[fortime] = forscore
    print "backward:"    
    for backtime in xrange(nexttime-1,pretime+1,-1):    
        assert backhist[backtime] == None and backhist[backtime+1] != None
        objclass = GreedyObjFunc(G,smodel,"backward",[backhist[backtime+1]])
        objmethod = getattr(objclass,"get{0}ObjFunc".format(algoinfo["objmethod"].capitalize()))
        restrict = getRestrictions(prestate,backhist[backtime+1],smodel)
        consts = genCoverCons(G,smodel,backhist[backtime+1])
        backhist[backtime],backscore = [deepcopy(item) for item in inferBackward(curstates[0],G,smodel,algoinfo,infermode,runfolder)]
        #backhist[backtime],backscore = [deepcopy(item) for item in detNonmonoSubCoverMaxSingle(objmethod,deepcopy(backhist[backtime+1]),consts,smodel,restrict)]
        backhistscore[backtime] = backscore
    maxscore,maxbipoint,maxbistate = None,None,None
    for bipoint in xrange(pretime+1,nexttime):
        assert forhist[bipoint-1] != None and backhist[bipoint+1] != None and forhistscore[bipoint-1] != None and backhistscore[bipoint+1] != None
        objclass = GreedyObjFunc(G,smodel,"bidirection",[forhist[bipoint-1],backhist[bipoint+1]])
        objmethod = getattr(objclass,"get{0}ObjFunc".format(algoinfo["objmethod"].capitalize()))
        print "bifind",bipoint,pretime,nexttime
        if len(backhist[bipoint+1][Trace.SUSCEPTIBLE].difference(forhist[bipoint-1][Trace.SUSCEPTIBLE])) != 0:
           continue
        restrict = getRestrictions(forhist[bipoint-1],backhist[bipoint+1],smodel)
        consts = genCoverCons(G,smodel,backhist[bipoint+1])
        tempbistate,biscore = [deepcopy(item) for item in detNonmonoSubCoverMaxSingle(objmethod,deepcopy(backhist[bipoint+1]),consts,smodel,restrict)] 
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
    
    
def inferBackward(curstate,G,smodel,algoinfo,infermode,runfolder):
    """infer backwards
    Args:
       curstate:
       G:
       smodel:
       algoinfo:
       infermode:
       runfolder:
    Returns:
       history:
    """
    seenstate = deepcopy(curstate)
    ensemble = algoinfo["ensemble"]
    enscount = 1
    if ensemble:
       enscount = algoinfo["enscount"]
    allsol = []
    allstates = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    for index in xrange(enscount):
        curstate = deepcopy(seenstate)
        sol,curcount= {}, 0
        while not HistoryUtil.isStopCond(sol,G,smodel,"Perfect","MinCut"):
           curstate = runMinCut(curstate,smodel,G,runfolder)
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
        print "inference before fixing"
        for time in sorted(sol.keys()):
            print time,len(sol[time][Trace.SUSCEPTIBLE]),len(sol[time][Trace.EXPOSED]),len(sol[time][Trace.INFECTED]),len(sol[time][Trace.RECOVERED]) 
        print "inference after fixing"        
        sol = HistoryUtil.lastFixSolution(sol,G,smodel,"None",infermode)
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
          
isTest = True
def runner(algoinfo,curstates,G,smodel,infermode,resultfile,runfolder,obstimes,TESTMODE=False):
    """iterative minimum cut relaxation for seir model
       Args:
         algoinfo:
         curstates: current state knowledge
         G: Graph
         smodel: spreading model
         infermode: infermode
         resultfile:
         runfolder:
         obstimes:
         TESTMODE: 
       Returns:
         sol: solution
    """
    assert smodel == "seir"
    globals()["isTest"] = TESTMODE
    if len(curstates[0][Trace.SUSCEPTIBLE]) == G.number_of_nodes():
       print "Data does not have any information!!All nodes are susceptible"
       exit(1)
    assert checkUtil.checkPreTrace(curstates,obstimes) and sorted(obstimes) == obstimes and len(obstimes) == len(curstates) and obstimes[0] == 0
    temphistory = inferBackward(curstates[0],G,smodel,algoinfo,infermode,runfolder,{})[0]
    if infermode == "Spreader":
       if algoinfo["ensemble"]: 
          enscount = algoinfo["enscount"]
          assert enscount == len(temphistory)     
          history = deepcopy(DisEns.getEnsembleSolution(temphistory,infermode,curstates,obstimes,G,smodel))
       else:   
          history = deepcopy(temphistory)
    elif infermode == "History": 
       assert not temphistory.has_key(0)     
       inittime = -1*min(temphistory.keys())
       history = {}
       for time in temphistory.keys():
           history[time+inittime] = deepcopy(temphistory[time])
       for index in xrange(1,len(obstimes)):
           pretime,nexttime = obstimes[index-1:index+1]
           prestate,nextstate = curstates[index-1:index+1]
           if nexttime-pretime == 1: 
              continue 
           print "iter: ",index,pretime,nexttime
           temphistory = inferBetweenDP(deepcopy(prestate),deepcopy(nextstate),pretime,nexttime,G,smodel,algoinfo)
           for time in temphistory.keys():
               assert not history.has_key(time+inittime)
               history[time+inittime] = deepcopy(temphistory[time])
       for tindex in xrange(len(obstimes)):
           time = obstimes[tindex]
           assert not history.has_key(time+inittime)
           history[time+inittime] = deepcopy(curstates[tindex]) 
       print len(history.keys())
       print "runnner insis"              
       for time in sorted(history.keys()):
           for state in history[time].keys():
               print time,state,len(history[time][state])
       print history.keys()
       print inittime
       print obstimes        
       history = InputOutput.history2Times(history,curstates,smodel)
    if globals()["isTest"]:
       checkUtil.checkSolution(infermode,"None",history,G,smodel,"notdetailed")
    InputOutput.writeHistoryResult(history,infermode,resultfile)  
    

    
if __name__ == "__main__":
    with gzip.open(sys.argv[1],"rb") as infile: 
        algoinfo = cPickle.load(infile)       
        curstates = cPickle.load(infile)
        G = cPickle.load(infile)
        smodel = cPickle.load(infile)
        infermode = cPickle.load(infile)
        resultfile = cPickle.load(infile)
        runfolder = cPickle.load(infile)
        obstimes = cPickle.load(infile)
    runner(algoinfo,curstates,G,smodel,infermode,resultfile,runfolder,obstimes)
