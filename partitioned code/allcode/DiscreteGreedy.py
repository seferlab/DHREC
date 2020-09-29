import networkx as nx
import numpy as np
import scipy as sp
import scipy.optimize
import random
import math
import sys
sys.path.append("./lib")
import os
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


class DiscreteGreedyTest():
     #Discrete Greedy Test Class
     def __init__():
         return

     @staticmethod
     def testForwardInference(forhist,prestate,nextstate):
         """tests forward reconstruction
         Args:
            forhist:
            prestate:
            nextstate:
         Returns:
            bool: true /false
         """ 
         sortedkeys = sorted(forhist.keys())
         for index in xrange(len(sortedkeys)-1):
             pretime,curtime = sortedkeys[index:index+2]  
             assert sum([len(forhist[pretime][state]) for state in forhist[pretime].keys()]) == sum([len(forhist[curtime][state]) for state in forhist[curtime].keys()]) 
         for time in forhist.keys():
             assert forhist[time] != None
         return True
     
     @staticmethod 
     def testVarMap(varmap,smodel,curstate):
         """tests varmap
         Args:
            varmap:
            smodel:
            curstate:
         Returns:
            bool: true/false
         """
         for node in varmap.keys():
             if varmap[node] == "i":
                assert node in curstate[Trace.RECOVERED] 
             elif  varmap[node] == "s":
                if smodel == "seir": 
                   assert node in curstate[Trace.EXPOSED]
                elif smodel in ["si","sir"]: 
                   assert node in curstate[Trace.INFECTED]   
             elif varmap[node] == "e":
                assert node in curstate[Trace.INFECTED] 
         for node in varmap.keys():
             if smodel == "si":
                assert varmap[node] in ["s"]
             elif smodel == "sir":
                assert varmap[node] in ["s","r"] 
             elif smodel == "seir":
                assert varmap[node] in ["s","e","r"]        
         return True

     @staticmethod
     def testAllSolution(history,G,smodel,curstate):    
         """test submodular solution of all
         Args:
            history:
            G:
            smodel:
            curstate:
         Returns:
            bool: true/false
         """
         assert checkUtil.checkCoverBetween(history[-1],curstate,smodel,G)
         histtimes = sorted(history.keys())   
         for tindex in xrange(1,len(histtimes)):
             pretime,curtime = histtimes[tindex-1:tindex+1]
             assert checkUtil.checkCoverBetween(history[pretime],history[curtime],smodel,G)
         return True    

     @staticmethod
     def testSolutionFeasible(X,item2cons,consts):
         """tests solution feasibility
         Args:
            X:
            item2cons
            consts: 
         Returns:
            bool: true/false
         """
         for item in X:
             assert isSolFeasible(item,item2cons,consts)
             if item2cons.has_key(item):
                 for cindex in item2cons[item]:
                     consts[cindex].removeItem(item)
         return True
                 
     @staticmethod 
     def testStepSolution(varmap,smodel,curstate,consts,item2cons,newstate,storeitem2cons,storeconsts,X):
         """test submodular solution of one step
         Args:
            varmap:
            smodel:
            curstate:
            consts:
            item2cons:
            newstate:
            storeitem2cons:
            storeconsts:
            X:
         Returns:
            bool: true/false
         """
         assert DiscreteGreedyTest.testVarMap(varmap,smodel,curstate)
         if smodel in ["si","sir"]:
            assert len(consts) == len(curstate[Trace.INFECTED])
         elif smodel == "seir":
            assert len(consts) == len(curstate[Trace.EXPOSED]) 
         assert len(item2cons.keys()) <= len([node for state in [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED] for node in curstate[state]])
         DiscreteGreedyTest.testSolutionFeasible(X,deepcopy(storeitem2cons),deepcopy(storeconsts))
         assert len([node for state in newstate.keys() for node in newstate[state]]) == len([node for state in curstate.keys() for node in curstate[state]])
         return True
     

class GreedyObjFunc():
     #Objective Function class
     def __init__(self,G,smodel,direction,inputstates):
         self.G = nx.DiGraph(G)
         self.smodel = smodel
         assert direction in ["backward","bidirection"]
         self.direction = direction
         self.states = HistoryUtil.getStates(smodel)
         self.statekey, self.exotype, self.transkeys = self.getTransInfo(smodel)
         if direction == "backward":
            self.curstate = deepcopy(inputstates[0]) #higher t
         elif direction == "bidirection":
            self.prestate = deepcopy(inputstates[0]) #lower t
            self.curstate = deepcopy(inputstates[1]) #higher t
         self.LOGZERO = -1000000000.0
         self.ZERO = 0
         self.TESTMODE = True

     def getTransInfo(self,smodel):
         """
         Args:
            smodel:
         Returns:
            statekey:
            transkeys:
         """
         if smodel in ["si", "sir", "sis"]:
            statekey = Trace.INFECTED
            exotype = "s2i"
         elif smodel == "seir":
            statekey = Trace.EXPOSED
            exotype = "s2e"
         if smodel == "seir":
            transkeys = ["s2s","s2e","e2e","e2i","i2i","i2r","r2r"] 
         elif self.smodel == "si":
            transkeys = ["s2s","s2i","i2i"]  
         elif smodel == "sir":
            transkeys = ["s2s","s2i","i2i","i2r","r2r"]  
         elif smodel == "sis":
            transkeys = ["s2s","s2i","i2s","i2i"]
         return statekey,exotype,transkeys
            
     def getUnitProb(self,info):
         """
         """
         return 1.0 - math.exp(-1.0*info[1][0])

     def getLogObjFunc(self,elemset):
         """get logarithm objective method
         Args:
            elemset: current elements
         Returns:
            score:
         """
         return self.getObjFunc(elemset,"log")
         
     def getNormalObjFunc(self,elemset):
         """get normal objective method
         Args:
            elemset: current elements
         Returns:
            score:
         """
         return self.getObjFunc(elemset,"normal")
     
     def getNodes(self,elemset):
         """gets useful nodes
         Args:
            elemset:
         Returns:
            nodes:
         """
         transnodes = {"s2e":set(),"e2e":set(),"s2i":set(),"i2i":set(),"s2s":set(),"i2s":set(),"e2i":set(),"i2r":set(),"r2r":set()}
         for key in self.transkeys:
             state1,state2 = key.split("2")
             transnodes[key] = set(elemset[state1].keys()).intersection(set(self.curstate[state2])) 
         for key in transnodes.keys():
             if len(transnodes[key]) == 0:
                del transnodes[key]
         if self.TESTMODE:
            nodeset = {"s":set(),"e":set(),"i":set(),"r":set()} 
            for key in transnodes: 
                state1,state2 = key.split("2")
                assert len(nodeset[state2].intersection(transnodes[key])) == 0 
                nodeset[state2] |= set(transnodes[key])
            for state in nodeset.keys():
                assert len(nodeset[state].symmetric_difference(self.curstate[state])) == 0    
         return transnodes,elemset        

     def getNodeTransScore(self,state1,state2,s2qnodes,elemset):
         """ gets transition scores for nodewise scores(i2s,e2e,e2i,i2r,r2r,i2i)
         Args:
            state1:
            state2: 
            s2qnodes: state s(state1) to state q(state2) nodes
            elemset:
         Returns:
            objval:   
         """
         if self.smodel == "si" and state1 == Trace.INFECTED and state2 == Trace.INFECTED: 
            return 0.0
         if state1 == Trace.RECOVERED and state2 == Trace.RECOVERED:
            return 0.0 
         objval = 0.0
         if state1 == state2:
            if state1 == Trace.EXPOSED: 
               transstr = "{0}2{1}".format(state1.capitalize(),Trace.INFECTED.capitalize())
            elif state1 == Trace.INFECTED:
               if self.smodel in ["seir","sir"]:
                  transstr = "{0}2{1}".format(state1.capitalize(),Trace.RECOVERED.capitalize())    
               elif self.smodel == "sis":
                  transstr = "{0}2{1}".format(state1.capitalize(),Trace.SUSCEPTIBLE.capitalize())     
         else:    
            transstr = "{0}2{1}".format(state1.capitalize(),state2.capitalize())
         for node in s2qnodes:
             if state1 != state2:
                localval = math.log(self.getUnitProb(self.G.node[node][getattr(Trace,transstr)]))
             else:
                localval = math.log(1.0 - self.getUnitProb(self.G.node[node][getattr(Trace,transstr)])) 
             objval += localval
         return objval
     
     def gets2sScore(self,s2snodes,elemset):
         """gets s to s transition
         Args:
            s2snodes:
            elemset:
         Returns:
            objval:
         """
         objval = 0.0 
         for node in s2snodes:
             curobjval = reduce(lambda x, y: x+y, [math.log(1.0 - self.getUnitProb(self.G[neighnode][node][self.exotype])) for neighnode in self.G.predecessors(node) if neighnode in elemset[Trace.INFECTED].keys()] + [0.0])
             objval += curobjval
         return objval
     
     def gets2eScore(self,s2inodes,s2enodes,elemset,mode):
         """gets s to e transition
         Args:
            s2inodes:
            s2enodes:
            elemset:
            mode:
         Returns:
            objval:
         """
         objval = 0.0
         for node in s2inodes.union(s2enodes):
             curobjval = 1.0 - reduce(lambda x, y: x*y, [1.0 - self.getUnitProb(self.G[neighnode][node][self.exotype]) for neighnode in self.G.predecessors(node) if neighnode in elemset[Trace.INFECTED].keys()] + [1.0])
             if curobjval == 0.0:
                if mode == "log":
                   return self.LOGZERO 
                elif mode == "normal":
                   return self.ZERO  
             else:
                objval += math.log(curobjval)
         return objval

     def getForwardScore(self,trans):
         """gets forward coefs
         Args:
           trans:
         Returns:
           objval:
         """
         objval = 0.0
         if self.smodel in ["si","sir"]:
            exotrans = "s2i"
         elif self.smodel == "seir":
            exotrans = "s2e"
         for node in trans.keys():
             prestate,nextstate = trans[node]
             if prestate == Trace.SUSCEPTIBLE and (nextstate == Trace.INFECTED or nextstate == Trace.EXPOSED):
                prob = 1.0 - reduce(lambda x, y: x*y, [1.0 - self.getUnitProb(self.G[neighnode][node][exotrans]) for neighnode in set(self.G.predecessors(node)).intersection(self.prestate[Trace.INFECTED])] + [1.0])
                if abs(prob-0.0) < 0.00000001:
                   objval += self.LOGZERO
                else:   
                   objval += math.log(prob)
             elif prestate == Trace.SUSCEPTIBLE and nextstate == Trace.SUSCEPTIBLE:
                prob = reduce(lambda x, y: x*y, [1.0 - self.getUnitProb(self.G[neighnode][node][exotrans]) for neighnode in set(self.G.predecessors(node)).intersection(self.prestate[Trace.INFECTED])] + [1.0])
                objval += math.log(prob)
             elif prestate == Trace.INFECTED and nextstate == Trace.INFECTED and self.smodel in ["sir","seir"]:
                objval += math.log(1.0-self.getUnitProb(self.G.node[node]["i2r"]))
             elif prestate == Trace.INFECTED and nextstate == Trace.RECOVERED:
                objval += math.log(self.getUnitProb(self.G.node[node]["i2r"]))
             elif prestate == Trace.EXPOSED and nextstate == Trace.INFECTED:
                objval += math.log(self.getUnitProb(self.G.node[node]["e2i"]))
             elif prestate == Trace.EXPOSED and nextstate == Trace.EXPOSED:
                objval += math.log(1.0-self.getUnitProb(self.G.node[node]["e2i"]))
         return objval
     
     def getObjFunc(self,elemset,mode):
         """returns MLE objective function in original form
         Args:
            elemset: (hash of hash) set of elements in the set
                 corr to susceptible for cur infected ones
                 corr to infected ones for cur recovered ones
            mode: log or normal
         Returns:
            objval: objective function value
         """
         totobjval = 0.0
         transnodes,newelemset = self.getNodes(elemset)
         objval = 0.0
         for key in transnodes.keys():
             state1,state2 = key.split("2")
             if key == "s2i":
                objval += self.gets2eScore(transnodes[key],set(),newelemset,mode)
             elif key == "s2e":
                objval += self.gets2eScore(set(),transnodes[key],newelemset,mode)
             elif key == "s2s":
                objval += self.gets2sScore(transnodes[key],newelemset)
             else:
                objval += self.getNodeTransScore(state1,state2,transnodes[key],newelemset)
         if mode == "log":    
            totobjval += objval
         elif mode == "normal":
            totobjval += math.exp(objval)
         if self.direction == "bidirection":
            premap = {}
            for state in self.prestate.keys():
                for node in self.prestate[state]:
                    premap[node] = state   
            trans = {}
            for state in elemset.keys():
                for node in elemset[state].keys():
                    if abs(elemset[state][node]-1.0) <= 0.00000001:
                       trans[node] = (premap[node],state) 
            assert len(set(trans.keys()).symmetric_difference(set(node for state in self.curstate.keys() for node in self.curstate[state]))) == 0 
            objval = self.getForwardScore(trans)
            if mode == "log":    
               totobjval += objval
            elif mode == "normal":
               totobjval += math.exp(objval)              
         return totobjval
         
def isSolFeasible(item,item2cons,consts):    
    """checks whether solution satisfies the constraints
      Args:
         item: list of elements 
         item2cons: elems to their covering constraints index mapping
         consts: constraints
      Returns:
         flag: boolean
    """
    if not item2cons.has_key(item):
       return True  
    for cindex in item2cons[item]:
        if not consts[cindex].isFeasible():
           return False  
    return True 

def getMap(curstate,smodel):
    """ returns variables to be added
    Args:
       curstate: current state
       smodel: spreading model
    Returns:
       varmap: node to variable map 
    """
    varmap = {}
    if smodel == "si":
       for node in curstate[Trace.INFECTED]:
           varmap[node] = Trace.SUSCEPTIBLE 
    elif smodel == "sis":
       for node in curstate[Trace.SUSCEPTIBLE]:
           varmap[node] = Trace.SUSCEPTIBLE
       for node in curstate[Trace.INFECTED]:
           varmap[node] = Trace.SUSCEPTIBLE     
    elif smodel == "seir":
       for node in curstate[Trace.EXPOSED]:
           varmap[node] = Trace.SUSCEPTIBLE
       for node in curstate[Trace.INFECTED]:
           varmap[node] = Trace.EXPOSED
       for node in curstate[Trace.RECOVERED]:
           varmap[node] = Trace.RECOVERED
    elif smodel == "sir":
       for node in curstate[Trace.INFECTED]:
           varmap[node] = Trace.SUSCEPTIBLE
       for node in curstate[Trace.RECOVERED]:
           varmap[node] = Trace.RECOVERED
    return varmap

def prepNewState(elemset,curstate,varmap,smodel):
    """makes state of curstate and current elements
    Args:
      elemset: current items 
      curstate: current state
      varmap:
      smodel: spreading model
    Returns:
      statemap: new state map after elemset modifications
    """
    statemap = deepcopy(curstate)
    if smodel in ["sir","seir"]:
       statemap[Trace.INFECTED] |= set(curstate[Trace.RECOVERED])
       statemap[Trace.RECOVERED] = set()
    if globals()["isTest"]:
       DiscreteGreedyTest.testVarMap(varmap,smodel,curstate)
    for node in elemset:
        if varmap[node] == Trace.SUSCEPTIBLE:
           if smodel in ["si","sir"]:
              assert node in curstate[Trace.INFECTED]
              statemap[Trace.INFECTED].remove(node)
           elif smodel == "seir":
              assert node in curstate[Trace.EXPOSED]
              statemap[Trace.EXPOSED].remove(node)
           elif smodel == "sis":
              assert node in curstate[Trace.INFECTED] or node in curstate[Trace.SUSCEPTIBLE]
              if node in curstate[Trace.INFECTED]:
                 statemap[Trace.INFECTED].remove(node)
           statemap[Trace.SUSCEPTIBLE].add(node)
        elif varmap[node] == Trace.EXPOSED:
           assert node in curstate[Trace.INFECTED]
           statemap[Trace.INFECTED].remove(node)
           statemap[Trace.EXPOSED].add(node)
        elif varmap[node] == Trace.RECOVERED:
           assert node in curstate[Trace.RECOVERED]
           statemap[Trace.RECOVERED].add(node)
           statemap[Trace.INFECTED].remove(node)
    return statemap           

def getConsInvMap(consts):
    """ returns inverse map 
    """
    node2cons = {} 
    for index in xrange(len(consts)):
        for node in consts[index].items:
            node2cons.setdefault(node,set()).add(index)
    return node2cons
       
def randNonmonoSubCoverMaxDouble(objmethod,curstate,consts,smodel):
    """Randomized Nonmonotone submodular maximization with additional covering constraints
      Args:
         objmethod: objective function to be evaluated
         curstate : current state of nodes
         consts: covering constraints list(list of constraint object)
         smodel: spreading model
      Returns:
         X: submodular maximizing elements
    """
    varmap = getMap(curstate,smodel)
    item2cons = getConsInvMap(consts)
    itemlist = varmap.keys()
    while True:
       random.shuffle(itemlist)
       X, Y = set(), set(itemlist)
       Xstate,Ystate = prepNewState(X,curstate,varmap,smodel), prepNewState(Y,curstate,varmap,smodel)
       Xstate = {state : {node: 1.0 for node in Xstate[state] } for state in Xstate.keys() }
       Ystate = {state : {node: 1.0 for node in Ystate[state] } for state in Ystate.keys() }
       prex,prey = objmethod(Xstate), objmethod(Ystate)
       for item in itemlist:
           if not isSolFeasible(item,item2cons,consts):
              Y.remove(item)
              prey = objmethod(Ystate)
              continue
           X.add(item)
           Y.remove(item)
           Xstate,Ystate = prepNewState(X,curstate,varmap,smodel), prepNewState(Y,curstate,varmap,smodel)
           Xstate = {state : {node: 1.0 for node in Xstate[state] } for state in Xstate.keys() }
           Ystate = {state : {node: 1.0 for node in Ystate[state] } for state in Ystate.keys() }
           curx,cury = objmethod(Xstate), objmethod(Ystate)
           ai = max(curx-prex, 0.0)
           bi = max(cury-prey, 0.0)
           flag = 0
           if ai == 0.0 and ai + bi == 0.0:
              if random.random() < 0.5:
                 flag = 1
              else:
                 flag = 2
           if flag == 1 or (flag ==0 and random.random() <= float(bi) / (ai+bi)):
              X.remove(item)
              prey = cury
           else:
              Y.add(item)
              prex = curx
              if item2cons.has_key(item):
                 for cindex in item2cons[item]:
                     consts[cindex].removeItem(item)
       assert len(X.difference(Y)) == 0
       break       
    newstate = prepNewState(X,curstate,varmap,smodel)
    return newstate


def getBestItem(mode,bestmode,X,itemlist,objmethod,curx,item2cons,consts,sideparams,remforbidX,addforbidX):
    """returns the best increasing item(either add or remove)
    Args:
       mode: add or remove
       bestmode: returns the best element or returns the one that is found
       X:
       itemlist:
       objmethod:
       curx: current objective
       item2cons:
       consts: 
       sideparams:
       remforbidX:
       addforbidX:
    Returns:
       bestitem:
       bestincr:
    """
    assert mode in ["add","remove"]
    assert bestmode in ["best","notbest"]
    [smodel,varmap,THRES,curstate]= sideparams
    bestitem, bestincr = None, 0  
    if mode == "add":
       tryset = list(set(itemlist).difference(addforbidX))
    elif mode == "remove":
       tryset = list(set(X).difference(remforbidX)) 
    if bestmode == "notbest":
       random.shuffle(tryset)
    for item in tryset:
        tempX = set(X)
        if mode == "add":
           tempX.add(item)
        elif mode == "remove":
           tempX.remove(item)
        for state in curstate.keys():
            assert type(curstate[state]) == set    
        tempXstate = prepNewState(tempX,curstate,varmap,smodel)
        tempXstate = {state : {node: 1.0 for node in tempXstate[state] } for state in tempXstate.keys() }
        tempinc = objmethod(tempXstate) - curx
        if tempinc < 0:
           continue 
        if curx == 0:
           continue  
        if (bestmode == "best" or (bestmode == "notbest" and abs(float(tempinc)/curx) >= THRES-1.0)) and isSolFeasible(item,item2cons,consts):   
           bestincr = tempinc
           bestitem = item
           if bestmode == "notbest":
              break 
    return bestitem,bestincr


def modifyVars(X,founditem,foundincr,mode,consts,item2cons,itemlist,curx,sideparams):
    """modify variables
    Args:
       X:
       founditem:
       foundincr:
       mode:
       consts:
       item2cons:
       itemlist:
       curx:
       sideparams:
    Returns:
       X,itemlist,consts,curx:
    """ 
    assert mode in ["add","remove"]  
    [smodel,varmap,THRES,curstate]= sideparams
    if mode == "add":    
       X.add(founditem)
       itemlist.remove(founditem)
    elif mode == "remove":
       X.remove(founditem)
       itemlist.add(founditem) 
    Xstate = prepNewState(X,curstate,varmap,smodel)
    Xstate = {state : {node: 1.0 for node in Xstate[state] } for state in Xstate.keys() }
    curx += foundincr 
    if mode == "add":
       if item2cons.has_key(founditem):
          for cindex in item2cons[founditem]:
              consts[cindex].removeItem(founditem)
    elif mode == "remove":
       if item2cons.has_key(founditem):
          for cindex in item2cons[founditem]:
              consts[cindex].addItem(founditem)          
    return X,itemlist,consts,curx


def prepareSol(X,curx,sideparams,startitem2cons,startconsts,objmethod):
    """prepares solution
    Args: 
       X:
       curx:
       sideparams:
       startitem2cons:
       startconsts:
       objmethod:
    Returns:
       newstate:   
    """
    [smodel,varmap,MINTHRES,curstate] = sideparams
    difX = set(varmap.keys()).difference(X)
    difXstate = prepNewState(difX,curstate,varmap,smodel)
    difXstate = {state : {node: 1.0 for node in difXstate[state] } for state in difXstate.keys()}
    flag = True 
    for item in difX:
        if not isSolFeasible(item,startitem2cons,startconsts):
           flag = False
           break
        if startitem2cons.has_key(item):
           for cindex in startitem2cons[item]:
               startconsts[cindex].removeItem(item)         
    if not flag:
       newstate = prepNewState(X,curstate,varmap,smodel)
    else:
       curdifx = objmethod(difXstate)
       if curx >= curdifx:
          newstate = prepNewState(X,curstate,varmap,smodel)
       else:
          newstate = prepNewState(difX,curstate,varmap,smodel)
    return newstate 


def TwoSnapshotValid(highstate,lowstate,smodel):
    """are two snapshots valid?
    Args:
       highstate: high t 
       lowstate: low t
       smodel:
    Returns:
       bool: true/false
    """
    allstates = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    for index1 in xrange(len(allstates)):
        state1 = allstates[index1]
        for index2 in xrange(index1+1,len(allstates)):
            state2 = allstates[index2]
            if state1 == Trace.EXPOSED and state2 == Trace.INFECTED:
               continue 
            print state1,state2,len(highstate[state1].intersection(lowstate[state2]))
            assert len(highstate[state1].intersection(lowstate[state2])) == 0
    return True

def getRestrictions(lowtstate,hightstate,smodel):
    """get restrictions given higher and lower bounds
    Args:
       hightstate: high t
       lowtstate: low t
       smodel:
    Returns:
       restrict:   
    """
    assert TwoSnapshotValid(hightstate,lowtstate,smodel)
    restrict = {}
    if smodel in ["si","sir","seir"]:
       assert sum([len(hightstate[state]) for state in hightstate.keys()]) == sum([len(lowtstate[state]) for state in lowtstate.keys()])
       assert len(hightstate[Trace.SUSCEPTIBLE].difference(lowtstate[Trace.SUSCEPTIBLE])) == 0 #may not be true for bi??
    for state in [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]:
        for node in hightstate[state].intersection(lowtstate[state]):
            restrict[node] = state
    return restrict      


def detNonmonoSubCoverMaxSingle(objmethod,curstate,consts,smodel,restrict={}):
    """deterministic nonmonotone submodular maximization with additional covering constraint(Implements the algorithm by Feige) 
      Args:
         objmethod: objective function to be evaluated
         curstate : current state of nodes
         consts: covering constraints list(list of constraint object)
         smodel:
         restrict: variable restrictions 
      Returns:
         newstate: submodular maximizing elements
         probscore: 
    """ 
    EPSILON = 1.0
    varmap = getMap(curstate,smodel)  
    item2cons = getConsInvMap(consts)
    startitem2cons,startconsts = deepcopy(item2cons),deepcopy(consts)
    remforbidX = set() #items that can not be removed
    addforbidX = set() #items that can not be added 
    for node in set(restrict.keys()).intersection(set(varmap.keys())):
        if varmap[node] == restrict[node]:
           remforbidX.add(node) 
           if item2cons.has_key(node):
              for cindex in item2cons[node]:
                  consts[cindex].removeItem(node)
        else:
           addforbidX.add(node) 
           #if item2cons.has_key(node):
           #   for cindex in item2cons[node]:
           #       consts[cindex].addItem(node)       
    itemlist = set(varmap.keys())        
    MINTHRES = 1.0 + (EPSILON / len(itemlist)**2)
    X = deepcopy(remforbidX)
    itemlist -= remforbidX
    Xstate = prepNewState(X,curstate,varmap,smodel)
    Xstate = {state : {node: 1.0 for node in Xstate[state]} for state in Xstate.keys() }
    curx = objmethod(Xstate)
    sideparams = [smodel,varmap,MINTHRES,curstate]
    for state in curstate.keys():
        assert type(curstate[state]) == set  
        
    #initial best element
    founditem,foundincr = getBestItem("add","best",X,itemlist,objmethod,curx,item2cons,consts,sideparams,remforbidX,addforbidX)
    if founditem == None or foundincr < 0:
       newstate = prepNewState(X,curstate,varmap,smodel) #?newsol or newstate
       tempnewstate = {state : {node: 1.0 for node in newstate[state] } for state in newstate.keys() }
       probscore = objmethod(tempnewstate)
       return newstate,probscore
    else:
       X,itemlist,consts,curx = modifyVars(X,founditem,foundincr,"add",consts,item2cons,itemlist,curx,sideparams)
       
    while True: #iterate over other elements 
       if globals()["isTest"]:
          assert len(set(X).intersection(set(itemlist))) == 0 and len(X) + len(itemlist) == len(varmap.keys()) 
       founditem,foundincr = getBestItem("add","notbest",X,itemlist,objmethod,curx,item2cons,consts,sideparams,remforbidX,addforbidX)
       if founditem != None:
          X,itemlist,consts,curx = modifyVars(X,founditem,foundincr,"add",consts,item2cons,itemlist,curx,sideparams)
          continue
       founditem,foundincr = getBestItem("remove","notbest",X,itemlist,objmethod,curx,item2cons,consts,sideparams,remforbidX,addforbidX)
       if founditem != None:
          X,itemlist,consts,curx = modifyVars(X,founditem,foundincr,"remove",consts,item2cons,itemlist,curx,sideparams)
          continue
       break
    newstate = prepareSol(X,curx,sideparams,deepcopy(startitem2cons),deepcopy(startconsts),objmethod)
    tempnewstate = {state : {node: 1.0 for node in newstate[state] } for state in newstate.keys() }
    probscore = objmethod(tempnewstate)
    if globals()["isTest"]:   
       DiscreteGreedyTest.testStepSolution(varmap,smodel,curstate,consts,item2cons,newstate,startitem2cons,startconsts,X)  
    return newstate,probscore
    

def genCoverCons(G,smodel,curstate):
    """generates covering constraints out of G and curstate
       Args:
         G: Graph
         smodel: spreading model
         curstate: current states
       Returns:
         consts: list of constraints of type Constraint
    """    
    consts = []
    affectstates = [Trace.INFECTED,Trace.RECOVERED]
    affectnodes = set(anode for astate in affectstates for anode in curstate[astate])
    if smodel in ["si","sir"]:
       nodes = curstate[Trace.INFECTED]
    elif smodel == "sis":
       nodes = [node for state in curstate.keys() for node in curstate[state].keys()]
    elif smodel == "seir":
       nodes = curstate[Trace.EXPOSED]
    for node in nodes:
        items = set([node]).union(set(G.predecessors(node)).intersection(affectnodes))
        right = len(items)-1
        cons = Constraint(items,right,"packing")
        consts.append(cons)
    return consts


def recursiveGreedySubCover(algoinfo,curstate,G,smodel,infermode,runfolder,typemethod):
    """infers history by greedy submodular algoritm
       Args:
         algoinfo: dict including algo parameters
         curstate: current states knowledge(can be more than one)
         G: Graph
         smodel: spreading model
         infermode: infer mode
         runfolder:
         typemethod: single/double
       Returns:
         history: solution
    """
    seenstate = deepcopy(curstate)
    assert typemethod in ["single","double"]
    itercount,ensemble = algoinfo["iter"],algoinfo["ensemble"]
    enscount = 1
    if ensemble:
       enscount = algoinfo["enscount"]
    allsol = []
    for index in xrange(enscount):
        curcount,history = 0, {} 
        curstate = deepcopy(seenstate)
        while (itercount == 'None' and not HistoryUtil.isStopCond(history,G,smodel,"Perfect",typemethod)) or (itercount != 'None' and curcount < itercount):
           print curcount,len(curstate[Trace.SUSCEPTIBLE]), len(curstate[Trace.EXPOSED]), len(curstate[Trace.INFECTED]), len(curstate[Trace.RECOVERED])
           objclass = GreedyObjFunc(G,smodel,"backward",[curstate])
           objmethod = getattr(objclass,"get{0}ObjFunc".format(algoinfo["objmethod"].capitalize()))
           consts = genCoverCons(G,smodel,curstate)
           if typemethod == "single":
              curstate = detNonmonoSubCoverMaxSingle(objmethod, curstate, consts,smodel,{})[0]
           elif typemethod == "double":
              curstate = randNonmonoSubCoverMaxDouble(objmethod, curstate, consts,smodel)
           for state in curstate.keys():
               assert type(curstate[state]) == set   
           history[-1 * (curcount + 1)] = deepcopy(curstate)
           curcount += 1
        if globals()["isTest"]:
           DiscreteGreedyTest.testAllSolution(history,G,smodel,seenstate)  
           exacttimes = InputOutput.history2Times(history,seenstate,smodel)  
           assert checkUtil.checkDiffusionOrder(exacttimes,smodel,G)
        history = HistoryUtil.FixSameParts(history,G)
        history = HistoryUtil.lastFixSolution(history,G,smodel,"None",infermode)
        print "index"              
        for time in sorted(history.keys()):
            print time,len(history[time][Trace.SUSCEPTIBLE]),len(history[time][Trace.EXPOSED]),len(history[time][Trace.INFECTED]),len(history[time][Trace.RECOVERED])
        allsol.append(deepcopy(history))  
    if globals()["isTest"]:
       for sol in allsol:
           checkUtil.checkPreConversion(sol,G,smodel)  
           assert checkUtil.basicCheck(sol,G)      
    if ensemble:
       return allsol 
    elif infermode == "Spreader":
       return allsol[0][min(allsol[0].keys())], min(allsol[0].keys())
    elif infermode == "History":
       return allsol[0]
     
def recursiveGreedySubCoverDouble(algoinfo,curstate,G,smodel,infermode,runfolder):
    """infers history by greedy submodular algoritm(double)
       Args:
         algoinfo: dict including algo parameters
         curstates: current states knowledge(can be more than one)
         G: Graphx
         smodel: spreading model
         infermode: infer mode
         runfolder:
       Returns:
         history: solution
    """ 
    return recursiveGreedySubCover(algoinfo,curstate,G,smodel,infermode,runfolder,"double")


def recursiveGreedySubCoverSingle(algoinfo,curstate,G,smodel,infermode,runfolder):
    """infers history by greedy submodular algoritm (single)
       Args:
         algoinfo: dict including algo parameters
         curstates: current states knowledge(can be more than one)
         G: Graph
         smodel: spreading model
         infermode: infer mode
         runfolder:
       Returns:
         history: solution
    """ 
    return recursiveGreedySubCover(algoinfo,curstate,G,smodel,infermode,runfolder,"single")


def getNextStates(smodel):
    """gets next states
    Args:
       smodel:
    Returns:
       state2state:
    """
    if smodel == "si":
       state2state = {Trace.SUSCEPTIBLE: Trace.INFECTED,Trace.INFECTED:Trace.INFECTED}
    elif smodel == "sir":
       state2state = {Trace.SUSCEPTIBLE: Trace.INFECTED,Trace.INFECTED: Trace.RECOVERED,Trace.RECOVERED:Trace.RECOVERED}
    elif smodel == "seir":
       state2state = {Trace.SUSCEPTIBLE: Trace.EXPOSED,Trace.EXPOSED:Trace.INFECTED,Trace.INFECTED:Trace.RECOVERED,Trace.RECOVERED:Trace.RECOVERED}
    return state2state   

def getUnitProb(info):
    """
    """
    return 1.0 - math.exp(-1.0*info[1][0])


def getForwardCoefs(nodeset,curstate,G,smodel):
    """
    Args:
       nodeset:
       curstate:
       G:
       snodel
    Returns:
       n2s2coef:
    """         
    if smodel in ["si","sir"]:
       exotrans = "s2i"
    elif smodel == "seir":
       exotrans = "s2e"
    n2s2coef = {}    
    state2next = getNextStates(smodel)  
    susnodes = nodeset.intersection(curstate[Trace.SUSCEPTIBLE])
    elsenodes = nodeset.difference(susnodes)   
    nextstate = state2next[Trace.SUSCEPTIBLE] 
    nodemap = {node: state for state in [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED] for node in curstate[state]}
    for node in susnodes:
        n2s2coef[node] = {}
        prob = 1.0 - reduce(lambda x, y: x*y, [1.0 - getUnitProb(G[neighnode][node][exotrans]) for neighnode in set(G.predecessors(node)).intersection(curstate[Trace.INFECTED])] + [1.0])
        n2s2coef[node][nextstate] = prob
        n2s2coef[node][Trace.SUSCEPTIBLE] = 1.0 - prob 
    for node in elsenodes:
        n2s2coef[node] = {}
        foundstate = nodemap[node]
        nextstate = state2next[foundstate] 
        if foundstate == nextstate:
           n2s2coef[node][foundstate] = 1.0 
           continue
        transkey = "{0}2{1}".format(foundstate.lower(),nextstate.lower())
        n2s2coef[node][nextstate] = getUnitProb(G.node[node][transkey])
        n2s2coef[node][foundstate] = 1.0 - n2s2coef[node][nextstate]
    if globals()["isTest"]:    
       for node in n2s2coef.keys():
           assert len(n2s2coef[node].keys()) <= 2 
           if len(n2s2coef[node].keys()) == 1:
              assert node in curstate[n2s2coef[node].keys()[0]]   
           else:
              first,second = n2s2coef[node].keys()
              assert state2next[first] == second or state2next[second] == first
    return n2s2coef
 

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
    inferstate = {Trace.SUSCEPTIBLE:set(),Trace.EXPOSED:set(),Trace.INFECTED:set(),Trace.RECOVERED:set()}
    inferscore = 0.0
    allnodes = set(node for state in curstate.keys() for node in curstate[state])
    state2next = getNextStates(smodel)
    endstate = [state for state in state2next.keys() if state == state2next[state]][0]
    endnodes = set(curstate[endstate])   
    feasnodes = allnodes.difference(endnodes)
    node2s2coef = getForwardCoefs(feasnodes,curstate,G,smodel)
    remainset = set()
    inferstate[endstate] |= set(curstate[endstate])
    for state in set(state2next.keys()).difference(set([endstate])):
        for node in curstate[state]:
            if restrict.has_key(node):
               inferstate[restrict[node]].add(node)
               prob = node2s2coef[node][restrict[node]]
               if abs(prob-0.0) < 0.00000001:
                  inferscore += -1000000000.0
               else:   
                  inferscore += math.log(prob)
            else:
               remainset.add(node)
    for node in remainset:
        selectstate = max(node2s2coef[node],key=node2s2coef[node].get)
        inferstate[selectstate].add(node)
        inferscore += math.log(node2s2coef[node][selectstate]) 
    if globals()["isTest"]:
       assert checkUtil.checkCoverBetween(curstate,inferstate,smodel,G)
       assert sum([len(inferstate[state]) for state in inferstate.keys()]) == sum([len(curstate[state]) for state in curstate.keys()])         
    return inferstate,inferscore

       
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
        backhist[backtime],backscore = [deepcopy(item) for item in detNonmonoSubCoverMaxSingle(objmethod,deepcopy(backhist[backtime+1]),consts,smodel,restrict)]
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


isTest = False
def runner(algo,algoinfo,infermode,curstates,G,smodel,resultfile,runfolder,obstimes,TESTMODE=True):
    """infers history by greedy algo
       Args:
          algo: algorithm
          algoinfo: algo parameters
          infermode: infer mode
          curstates: given states information(must be 1 in greedy case)
          G: Graph
          smodel:
          resultfile: result file
          runfolder:
          obstimes:
          TESTMODE: test
       Returns:
          sol: 
    """
    globals()["isTest"] = TESTMODE
    assert checkUtil.checkPreTrace(curstates,obstimes) and sorted(obstimes) == obstimes and len(obstimes) == len(curstates) and obstimes[0] == 0
    allstates = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    methodname = "recursive{0}".format(algo)
    temphistory = globals()[methodname](algoinfo,deepcopy(curstates[0]),G,smodel,infermode,runfolder)
    #for time in temphistory.keys():
    #    print time,len(temphistory[time][Trace.SUSCEPTIBLE]), len(temphistory[time][Trace.EXPOSED]), len(temphistory[time][Trace.INFECTED]), len(temphistory[time][Trace.RECOVERED])
    if infermode == "Spreader":
       if algoinfo["ensemble"]: 
          enscount = algoinfo["enscount"]
          assert enscount == len(temphistory)     
          history = deepcopy(DisEns.getEnsembleSolution(temphistory,infermode,curstates,obstimes,G,smodel))
          print "write history:"
          print history
       else:   
          history = deepcopy(temphistory)
    elif infermode == "History": 
       if False: 
        print obstimes 
        print temphistory.keys()
        print "curstae"
        for state in curstates[0].keys():
           print state,len(curstates[0][state]) 
        print "temphist"
        for state in temphistory[0].keys():
           print state,len(temphistory[0][state])     
        for state in curstates[0].keys():
           assert len(curstates[0][state]) == len(temphistory[0][state])
       assert not temphistory.has_key(0)     
       #del temphistory[0]  
       #temphistory = {0:set()} #? 
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
       assert len(set(history.keys()).symmetric_difference(range(0,max(history.keys())+1))) == 0    
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
        algo = cPickle.load(infile)
        algoinfo = cPickle.load(infile)
        infermode = cPickle.load(infile)
        curstates = cPickle.load(infile)
        G = cPickle.load(infile)
        smodel = cPickle.load(infile)
        resultfile = cPickle.load(infile)
        runfolder = cPickle.load(infile)
        obstimes = cPickle.load(infile)
    runner(algo,algoinfo,infermode,curstates,G,smodel,resultfile,runfolder)
