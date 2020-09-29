import networkx as nx
import numpy as np
import scipy as sp
import scipy.optimize
import random
import time
import math
import sys
sys.path.append("./lib")
import os
import itertools
from copy import deepcopy
import myutilities as myutil
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from HistoryEnsemble import DiscreteNoneEnsemble as DisEns
from DiscreteDistribution import DiscreteDistribution as Dist
import HistoryUtil
import checkUtil
from Qsapmin import QsapMinTest
import boundInferUtil


class QsapMaxTest():
     #Qsap Max Test Class
     def __init__():
         return

     @staticmethod
     def testReturnValues(retvalues,node2vars):
         """test return values
         """
         singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
         doublevars = [var for var in retvalues.keys() if var.find("_") != -1]
         for node in node2vars.keys():
             mysum = sum([retvalues[var] for var in node2vars[node] if retvalues.has_key(var)])
             assert abs(mysum - 1.0) <= 0.00001
         fracvals = set([retvalues[singlevar] for singlevar in singlevars]).difference(set([0.0,0.5,1.0]))
         assert len(fracvals) <= 2
         
     @staticmethod
     def testTransformCoefs(modlinear,modquad,node2vars):
         """tests transformed coefs
         Args:
            modlinear: 
            modquad: 
            node2vars: 
         Returns:
            bool: true/false 
         """
         seenvars = [var for node in node2vars.keys() for var in node2vars[node]]
         for val in modquad.values():
             assert val >= 0.0
         for val in modlinear.values():
             assert val >= 0.0  
         for var1,var2 in modquad.keys():
             assert not modquad.has_key((var2,var1)) 
             assert modquad[(var1,var2)] > 0.0 and var1 in seenvars and var2 in seenvars 
             node1 = int(var1.split("?")[0].replace("x",""))
             node2 = int(var2.split("?")[0].replace("x","")) 
             assert node2vars.has_key(node1) and node2vars.has_key(node2)
         for node in node2vars.keys():
             assert len(node2vars[node]) > 1
         for var in modlinear.keys():
             assert var in seenvars   
         for var in modlinear.keys():
             node = int(var.split("?")[0].replace("x",""))
             assert node2vars.has_key(node)      
         return True
     
     @staticmethod
     def testPreRun(algo,linear,quad):
         """prerun test
         Args:
           algo:
           linear:
           quad:
         Returns:
           bool: true/false
         """    
         if algo == "pipage":   
            for val in linear.values():
                assert val >= 0.0
         elif algo == "greedy": 
            pass  
         return True
     
     @staticmethod
     def testPipage(retvalues,node2vars,binsol,linear,quad):
         """test pipage rounding
         """
         retval1 = estimateLpObjVal(linear,quad,node2vars,binsol,"binary")
         fracsol = {}
         for node in node2vars.keys():
             fracsol[node] = {}
             for var in node2vars[node]:
                 val = 0.0
                 if retvalues.has_key(var):
                    val = retvalues[var]
                 fracsol[node][var] = val
         retval2 = estimateLpObjVal(linear,quad,node2vars,fracsol,"frac")
         assert retval1 >= retval2

     @staticmethod    
     def testBasics(binsol,history,node2slots,node2vars,retvalues):
         """tests basic variables
         """    
         states = [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
         histnodes = set(node for time in history.keys() for state in states for node in history[time][state])
         solnodes = set(binsol.keys())
         assert len(histnodes.symmetric_difference(solnodes)) == 0
         assert len(solnodes.difference(set(node2slots.keys()))) == 0 and len(set(node2slots.keys()).difference(solnodes)) == 0  
         assert len(set(node2slots.keys()).intersection(set(node2vars.keys()))) == len(set(node2slots.keys())) and len(set(node2slots.keys()).intersection(solnodes)) == len(solnodes)  
         singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
         doublevars = [var for var in retvalues.keys() if var.find("_") != -1]
         onenodes = [int(var.split("?")[0].replace("x","")) for var in singlevars if retvalues[var] == 1.0]
          
         startnodes = set()
         for node in node2vars.keys():    
             for varname in node2vars[node]:
                 temptimes = HistoryUtil.var2Times(varname)
                 if temptimes.has_key(Trace.INFECTED) and temptimes[Trace.INFECTED] == 0:
                    startnodes.add(node)                 
         difset = set(history[0][Trace.INFECTED]).difference(startnodes)  
         assert len(difset) == 0
         #for item in difset:
         #    temptimes = HistoryUtil.var2Times(node2vars[node])
         #    print item,temptimes,varname      
         
     @staticmethod
     def testLpSolution(paraminfo,solinfo):
         """lp solution test
         Args:
            paraminfo:
            solinfo:
         Returns:
            bool: true/false
         """
         [node2slots,node2vars,modcurstates,quad,linear,modquad,modlinear,G,smodel] = paraminfo
         [retvalues,binsol,history] = solinfo
         states = [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
         QsapMaxTest.testBasics(binsol,history,node2slots,node2vars,retvalues)
         QsapMaxTest.testTransformCoefs(modlinear,modquad,node2vars)
         QsapMaxTest.testReturnValues(retvalues,node2vars)
         QsapMaxTest.testPipage(retvalues,node2vars,binsol,modlinear,modquad)
         return True
     
     @staticmethod
     def testGreedyModifiedCoefs(linear,quad,modlinear,modquad):
         """checks modified coefs of greedy algo
         Args:
           linear:
           quad:
           modlinear:
           modquad:
         Returns:
           bool: true/false
         """
         for var1,var2 in quad.keys():
             assert not quad.has_key((var2,var1))
         for var1,var2 in modquad.keys():
             assert not modquad.has_key((var2,var1))
         neglinear = deepcopy(modlinear)    
         for var1,var2 in modquad.keys():
             if not neglinear.has_key(var1):
                neglinear[var1] = 0.0
             neglinear[var1] -= modquad[(var1,var2)]
             if not neglinear.has_key(var2):
                neglinear[var2] = 0.0
             neglinear[var2] -= modquad[(var1,var2)]
         #print "min after modification",min(neglinear.values())
         #print "sum info"
         #print sum(quad.values())
         #print sum(linear.values())
         #print sum(modlinear.values())
         #print abs(sum(quad.values()) + sum(linear.values()) - sum(modlinear.values()))
         assert abs(sum(quad.values()) + sum(linear.values()) - sum(modlinear.values())) <= 0.001 
         return True
     
     @staticmethod
     def testGreedyResults(linear,quad,modlinear,modquad,cursol,optval,node2vars):
         """check results of greedy algo by comparing etc
         Args:
            linear:
            quad:
            modlinear:
            modquad:
            cursol:
            optval:
            node2vars:
         Returns:
           bool: true/false 
         """
         solval = sum([linear[var] for node,var in cursol if linear.has_key(var)]) 
         solval += sum([quad[(var1,var2)] for node1,var1 in cursol for node2,var2 in cursol if quad.has_key((var1,var2))])
         solval2 = sum([modlinear[var] for node,var in cursol if modlinear.has_key(var)]) 
         solval2 -= sum([modquad[(var1,var2)] for node1,var1 in cursol for node2,var2 in cursol if modquad.has_key((var1,var2))])
         assert abs(optval - solval) <= 0.000001 and abs(optval - solval2) <= 0.000001 and abs(solval - solval2) <= 0.000001
         assert len(cursol) == len(node2vars.keys()) 
         for node,var in cursol.keys():
             assert node2vars.has_key(node) and var in node2vars[node]
         return True


def estimateGreedyModCoefs(linear,quad,node2vars):
    """estimates modified coefs for greedy algo
    Args:
       linear:
       quad:
       node2vars:
    Returns:
       modlinear:
       modquad:
    """
    modquad = {}
    modlinear = deepcopy(linear)
    for nodevar1,nodevar2 in quad.keys():
        modlinear.setdefault(nodevar1,0.0)
        modlinear[nodevar1] += quad[(nodevar1,nodevar2)]
        node2 = int(nodevar2.split("?")[0].replace("x",""))
        for remain in node2vars[node2] - set([nodevar2]):
            if modquad.has_key((nodevar1,remain)):
               quadkey = (nodevar1,remain)
            elif modquad.has_key((remain,nodevar1)):
               quadkey = (remain,nodevar1)
            else:
               quadkey = (nodevar1,remain) 
               modquad[quadkey] = 0.0
            modquad[quadkey] += quad[(nodevar1,nodevar2)]
    return modlinear,modquad     


def greedyMaxGroupCoverage(node2vars,linear,quad,maxtime):
    """implements greedy max group coverage algorithm for max case (1/4 approximation) 
    Args:
       node2vars:
       linear:
       quad:
       maxtime:
    Returns:
       history:
    """
    if globals()["isTest"]:
       assert QsapMaxTest.testPreRun("greedy",linear,quad)
    
    modlinear,modquad = estimateGreedyModCoefs(linear,quad,node2vars)
    cursol = set((node,list(node2vars[node])[0]) for node in node2vars.keys() if len(node2vars[node]) == 1)
    curnodes = [node for node in node2vars.keys() if len(node2vars[node]) > 1]
    random.shuffle(curnodes)
    optval = 0.0
    if len(cursol) != 0:
       optval = sum([modlinear[varname] for node,varname in cursol if modlinear.has_key(varname)])
       optval -= sum([modquad[(var1,var2)] for node1,var1 in cursol for node2,var2 in cursol if modquad.has_key((var1,var2))])
       
    for node in curnodes:
        maxnodevar,maxinc = None, -1000000000.0
        for varname in node2vars[node]:
            inc = 0.0
            if modlinear.has_key(varname):
               inc = modlinear[varname]
            for node1,var1 in cursol:
                if modquad.has_key((varname,var1)):
                   inc -= modquad[(varname,var1)]
                elif modquad.has_key((var1,varname)):
                   inc -= modquad[(var1,varname)] 
            if inc >= maxinc:
               maxnodevar = (node,varname)
               maxinc = inc  
        assert maxnodevar != None      
        cursol.add(maxnodevar)
        optval += maxinc
        
    while len(curnodes) != 0:
        maxnodevar,maxinc = None, -1000000000.0
        for node in curnodes:
            for varname in node2vars[node]:
                inc = 0.0
                if modlinear.has_key(varname):
                   inc = modlinear[varname]
                for node1,var1 in cursol:
                    if modquad.has_key((varname,var1)):
                       inc -= modquad[(varname,var1)]
                    elif modquad.has_key((var1,varname)):
                       inc -= modquad[(var1,varname)] 
                if inc >= maxinc:
                   maxnodevar = (node,varname)
                   maxinc = inc  
        assert maxnodevar != None      
        curnodes.remove(maxnodevar[0])
        cursol.add(maxnodevar)
        optval += maxinc  
    assert optval >= 0.0
    history = {time: {state: set() for state in [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]} for time in xrange(maxtime+1)}
    for node,var in cursol:
        times = HistoryUtil.var2Times(var)
        for state in times.keys():
            history[times[state]][state].add(node)
    if globals()["isTest"]:
       QsapMaxTest.testGreedyResults(linear,quad,modlinear,modquad,cursol,optval,node2vars)
       QsapMaxTest.testGreedyModifiedCoefs(linear,quad,modlinear,modquad)
    return history


def preprocessLp(linear,quad,node2vars):
    """preprocess for lp
    Args:
       linear:
       quad:
       node2vars:
    Returns:
       modlinear:
       modquad: 
    """
    presolvars = set(list(node2vars[node])[0] for node in node2vars.keys() if len(node2vars[node]) == 1 )
    modlinear, modquad = deepcopy(linear),deepcopy(quad)
    if len(presolvars) != 0:
       for var1,var2 in modquad.keys():
           if var1 in presolvars and var2 in presolvars:
              del modquad[(var1,var2)]
           elif var1 in presolvars:
              modlinear.setdefault(var2,0.0)
              modlinear[var2] += quad[(var1,var2)]
              del modquad[(var1,var2)]
           elif var2 in presolvars:
              modlinear.setdefault(var1,0.0)
              modlinear[var1] += quad[(var1,var2)]
              del modquad[(var1,var2)]
       for var in presolvars:
           node = int(var.split("?")[0].replace("x",""))
           del node2vars[node]
       for var in modlinear.keys():
           node = int(var.split("?")[0].replace("x",""))
           if not node2vars.has_key(node):
              del modlinear[var]
    return modlinear,modquad
             
        
def genMaxBoundConst(linear,quad,node2vars):
    """generates Maximization problem bound constraints
    Args:
       linear:
       quad:
       node2vars:
    Returns:
       boundstr:   
    """
    boundstr = "Bounds\n"
    boundstr += "\n".join(["0 <= {0} <= 1.0".format(varname) for node in node2vars.keys() for varname in node2vars[node]]) + "\n"
    boundstr += "\n".join(["0 <= {0}_{1} <= 1.0".format(var1,var2) for var1,var2 in quad.keys()]) + "\n"
    return boundstr


def genMaxObj(linear,quad):
    """generates max problem objective function 
    Args:
      linear:
      quad:
    Returns:
      objstr:
    """   
    isPlus = lambda x: "+" if x >= 0 else " "
    objstr = "Maximize\n obj: "
    objstr += " ".join([" {0} %.10f {1} ".format(isPlus(linear[var1]),var1) %linear[var1]  for var1 in linear.keys()])       
    objstr += " ".join([" {0} %.10f {1} ".format(isPlus(quad[(var1,var2)]),var1+"_"+var2) %quad[(var1,var2)] for var1,var2 in quad.keys()])
    return objstr
       
def genMaxMainConst(node2vars,linear,quad,smodel,G):
    """generates max problem constraints 
    Args:
       node2vars: 
       linear:
       smodel:
       G:
       quad:
    Returns:
       constr:
    """     
    consstr = "Subject To \n"
    consstr += "\n".join([" + ".join([" 1.0 {0} ".format(varname) for varname in node2vars[node]]) + " = 1.0 " for node in node2vars.keys()]) + "\n"
    for var1,var2 in quad.keys():
        consstr += "{0}_{1} - {0} <= 0.0\n".format(var1,var2)
        node2 = int(var2.split("?")[0].replace("x",""))
        for tvar in node2vars[node]:
            consstr += "{0}_{1} + {2} <= 1.0\n".format(var1,var2,tvar)
    if False:
       #generates covering constraints
       if smodel in ["si","sir"]:
          sstate,rstate = Trace.INFECTED,Trace.INFECTED
       elif smodel == "seir":
          sstate,rstate = Trace.INFECTED,Trace.EXPOSED       
       for node in node2vars.keys():
           for var1 in node2vars[node]:
               times1 = HistoryUtil.var2Times(var1)
               if times1[Trace.INFECTED] > 0: 
                  arrs = [" 1.0 {0} ".format(var2) for neighnode in set(G.predecessors(node)).intersection(set(node2vars.keys())) for var2 in node2vars[neighnode] if HistoryUtil.checkAffect(HistoryUtil.var2Times(var2),times1,smodel,sstate,rstate)]
               consstr += " + ".join(arrs) + " -1.0 {0} >= 0.0\n".format(var1)

    return consstr 


def pipageRound(fracsol,node2vars,linear,quad):
    """pipage rounding
    Args:
       fracsol:
       node2vars:
       linear:
       quad:
    Returns:
       retsol:
    """
    for node in fracsol.keys():
        assert abs(sum(fracsol[node].values())-1) < 0.0001 
    curobj = estimateLpObjVal(linear,quad,node2vars,fracsol,"frac")
    nodefracmap = {node:[var for var in fracsol[node].keys() if fracsol[node][var] >= 0.0001 or abs(fracsol[node][var]-1) >= 0.00001] for node in fracsol.keys()}  
    for onenode in [node for node in fracsol.keys() for var in fracsol[node] if abs(fracsol[node][var] - 1) < 0.0001]:
        del nodefracmap[onenode]
    presum = sum([len(nodefracmap[node]) for node in nodefracmap.keys()])  
    count = 0
    while True:
       print "curiter {0} {1}".format(len(nodefracmap.keys()),curobj)
       for node in fracsol.keys():
           assert abs(sum(fracsol[node].values())-1) < 0.0001 
       if len(nodefracmap.keys()) == 0:
          break 
       count += 1
       node = random.choice(nodefracmap.keys())
       curvars = nodefracmap[node]
       assert len(curvars) >= 2
       random.shuffle(curvars)
       var1,var2 = curvars[0:2]
       eps1 = min(1.0-fracsol[node][var1],fracsol[node][var2])
       eps2 = -1*min(fracsol[node][var1],1.0-fracsol[node][var2])
       if random.random() < 0.5:
          epslist  = [0,eps1,eps2]
       else:
          epslist  = [0,eps2,eps1]
       flag = False
       for index in xrange(1,len(epslist)):
           preeps,cureps = epslist[index-1:index+1]
           fracsol[node][var1] += cureps - preeps
           fracsol[node][var2] += preeps - cureps
           tobj = estimateLpObjVal(linear,quad,node2vars,fracsol,"frac")
           if (tobj >= curobj and index < len(epslist)-1) or ((tobj >= curobj or abs(tobj-curobj) <= 0.00001) and index == len(epslist)-1):
              curobj = tobj
              flag = True
              break
       assert flag
       if abs(fracsol[node][var1]) <= 0.0001:
          nodefracmap[node].remove(var1) 
       elif abs(fracsol[node][var2]) <= 0.0001:
          nodefracmap[node].remove(var2)    
       if abs(fracsol[node][var1]-1) <= 0.0001 or abs(fracsol[node][var2]-1) <= 0.0001:
          del nodefracmap[node] 
    assert count <= presum  
    retsol = {}   
    for node in fracsol.keys():
        curvar = [var for var in fracsol[node] if abs(fracsol[node][var]-1) < 0.0001]  
        assert len(curvar) == 1
        retsol[node] = curvar[0]
    assert len(retsol) == len(node2vars.keys())    
    return retsol
   

def estimateLpObjVal(linear,quad,node2vars,sol,mode):
    """estimates objective value frac
    Args:
       linear:
       quad:
       node2vars:
       sol:
       mode: frac or binary
    Returns:
       objval:  
    """
    assert mode in ["frac","binary"]
    if mode == "binary":
       objval = sum([linear[var] for var in sol.values() if linear.has_key(var)])
       objval += sum([quad[(var1,var2)] for var1 in sol.values() for var2 in sol.values() if quad.has_key((var1,var2))])
    elif mode == "frac":
       objval = sum([linear[var]*sol[node][var] for node in sol.keys() for var in sol[node].keys() if linear.has_key(var)])
       objval += sum([quad[(var1,var2)]*sol[node1][var1]*sol[node2][var2] for node1 in sol.keys() for var1 in sol[node1].keys() for node2 in sol.keys() for var2 in sol[node2].keys() if quad.has_key((var1,var2))])
    return objval

    
def getVcounts(fracsol,qs):
    """gets v counts
    Args:
       fracsol:
       qs:
    Returns:
       V1:
       v2:
    """  
    qscount = {q:0 for q in qs}
    for node in fracsol.keys():
        for var in fracsol[node].keys():
            if fracsol[node][var] in qs:
               qscount[fracsol[node][var]] += 1
    return qscount[min(qs)],qscount[max(qs)]   


def randomRound(fracsol):
    """assign var to each node by randomized rounding
    Args:
       fracsol:
    Returns:
       fracsol:
    """
    for node in fracsol.keys():
        curvars = list(fracsol[node])
        random.shuffle(curvars)
        foundvar = None
        curval = 0.0
        p = random.random() 
        for var in curvars:
            curval += fracsol[node][var]
            if curval >= p:
               foundvar = var
               break
        if foundvar == None:
           foundvar = curvars[-1]      
        fracsol[node][foundvar] = 1.0
        fracsol[node] = {curvar: 0.0 for curvar in set(curvars).difference(set([foundvar]))}    
    return fracsol    


def twoStepRound(retvalues,node2vars,linear,quad):
    """two step rounding procedure based on two step pipage rounding
    Args:
       retvalues:
       node2vars:
       linear:
       quad:
    Returns:
       sol:
    """    
    fracsol = {}
    vals = set()
    for node in node2vars.keys():
        fracsol[node] = {}
        for var in node2vars[node]:
            val = 0.0
            if retvalues.has_key(var):
               val = retvalues[var]
            fracsol[node][var] = val
            vals.add(val)
    fracsol2 = deepcopy(fracsol)    
    if globals()["isTest"]: 
       for node in fracsol2.keys():
           assert abs(sum(fracsol2[node].values())-1) < 0.0001     
    sol1 = pipageRound(deepcopy(fracsol),node2vars,linear,quad)
    deltas = list(vals.difference(set([0.0,0.5,1.0])))
    if len(deltas) == 0:
       return sol1
    if len(deltas) == 1:
       deltas.append(1-deltas[0])     
    assert len(deltas) >= 2
    delta,revdelta = min(deltas),max(deltas)
    V1 = len([1 for node in fracsol2.keys() for val in fracsol2[node].values() if abs(val-delta) < 0.0001])
    V2 = len([1 for node in fracsol2.keys() for val in fracsol2[node].values() if abs(val-revdelta) < 0.0001])
    print "deltas"
    print delta
    print revdelta
    print "vinfo:"
    print V1
    print V2
    for node in fracsol2.keys():
        for var in fracsol2[node].keys():
            if abs(fracsol[node][var] - delta) < 0.0001:
               fracsol2[node][var] = min(1.0,delta + ((revdelta*V2)/float(V1)))
            elif abs(fracsol[node][var] - revdelta) < 0.0001:
               fracsol2[node][var] = max(0.0,revdelta*(1-(V1/float(V2))))           
    if globals()["isTest"]: 
       for node in fracsol2.keys():
           print abs(sum(fracsol2[node].values())-1) 
           assert abs(sum(fracsol2[node].values())-1) < 0.0001 
    sol2 = pipageRound(deepcopy(fracsol2),node2vars,linear,quad) 
    if globals()["isTest"]:
       assert len(sol1.keys()) == len(sol2.keys()) 
    if estimateLpObjVal(linear,quad,node2vars,sol1,"binary") >= estimateLpObjVal(linear,quad,node2vars,sol2,"binary"):
       return sol1
    else:
       return sol2  
        
 
def transformCoefs(linear,quad,node2vars):
    """transforms coefficients
    Args:
       linear:
       quad:
       node2vars:
    Returns:
       linear:
       quad:
    """
    posquadvals = [(var1,var2) for var1,var2 in quad.keys() if quad[(var1,var2)] > 0.0]
    negquadvals = [(var1,var2) for var1,var2 in quad.keys() if quad[(var1,var2)] < 0.0]
    newquad = {(var1,var2): quad[(var1,var2)] for var1,var2 in negquadvals}
    for var1,var2 in posquadvals:
        node2 = int(var2.split("?")[0].replace("x","")) 
        linear.setdefault(var1,0.0)
        linear[var1] += quad[(var1,var2)]
        for remvar in set(node2vars[node2]).difference(set([var2])):
            if newquad.has_key((remvar,var1)):
               curkey = (remvar,var1)
            else:  
               curkey = (var1,remvar)  
               newquad.setdefault(curkey,0.0)
            newquad[curkey] += -1.0*quad[(var1,var2)] 
    quad = {(var1,var2):-1*newquad[(var1,var2)] for var1,var2 in newquad.keys()}
    poslinear = {}    
    for var in linear.keys():
        if linear[var] < 0.0:
           poslinear[var] = -1*linear[var]
           continue 
        node = int(var.split("?")[0].replace("x","")) 
        for sidevar in set(node2vars[node]).difference(set([var])):
            poslinear.setdefault(sidevar,0.0)
            poslinear[sidevar] += linear[var]
    linear = deepcopy(poslinear)
    return linear,quad            


def convert2Sol(binsol,maxobstime):
    """converts retvalues to sol format(history format)
    Args:
       binsol:
       maxobstime:
    Returns:
       history:
    """
    states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    history = {time: {state:set() for state in states} for time in xrange(maxobstime+1)}
    for node in binsol.keys():
        varname = binsol[node]
        node = int(varname.replace("x","").split("?")[0])
        nodetimes = HistoryUtil.var2Times(varname)
        for state in nodetimes.keys():
            history[nodetimes[state]][state].add(node)
    return history


def runLp(modcurstates,smodel,G,runfolder,node2slots,node2vars,linear,quad):
    """runs lp for max problem
    Args:
       modcurstates:
       smodel:
       G:
       runfolder:
       node2slots:
       node2vars:
       linear:
       quad:
    Returns:
       history:
    """  
    if globals()["isTest"]:
       assert QsapMaxTest.testPreRun("pipage",linear,quad)
    modlinear, modquad = preprocessLp(linear,quad,node2vars)
    modlinear,modquad = transformCoefs(modlinear,modquad,node2vars)
    maxtime = max(modcurstates.keys())
    objstr = genMaxObj(modlinear,modquad)
    consstr = genMaxMainConst(node2vars,modlinear,modquad,smodel,G)
    boundstr = genMaxBoundConst(modlinear,modquad,node2vars)
    retvalues,objval,runtime = HistoryUtil.runHistoryCode(consstr,objstr,boundstr,runfolder,["x"],"qsapmax-lp","removeZeros")
    binsol = twoStepRound(retvalues,node2vars,modlinear,modquad)
    history = convert2Sol(binsol,maxtime)
    if globals()["isTest"]:
       paraminfo = [node2slots,node2vars,modcurstates,quad,linear,modquad,modlinear,G,smodel]
       solinfo = [retvalues,binsol,history] 
       QsapMaxTest.testLpSolution(paraminfo,solinfo)
    return history 

               
def paramCheck(algoinfo,G,smodel):
    """check parameters of algo
    Args:
       algoinfo:
       G:
       smodel:
    Returns:
       bool: true/false
    """
    transmap = {"si":Trace.S2I,"sir":Trace.S2I,"seir":Trace.S2E}
    if algoinfo["method"] == "greedy":
       for node1,node2 in G.edges():
           if G[node1][node2][transmap[smodel]][0] not in ["expo","geometric"]:
              return False
    return True

isTest = False
def runner(algoinfo,curstates,obstimes,G,smodel,infermode,resultfile,runfolder,TESTMODE=False):
    """uniform metric labeling relaxation for known interval
       Args:
         algoinfo:
         curstates: current state knowledge
         obstimes:
         G: Graph
         smodel: spreading model
         infermode: infermode
         resultfile:
         runfolder:
         TESTMODE:
       Returns:
         sol: solution
    """
    APPROX = "taylor"
    globals()["isTest"] = TESTMODE
    method = algoinfo["method"]
    assert method in ["greedy","pipage"]
    if not paramCheck(algoinfo,G,smodel):
       "Parameters and method does not match!!"
       exit(1) 
    assert checkUtil.checkPreTrace(curstates,obstimes)
    maxtime = max(obstimes)
    ensemble = algoinfo["ensemble"]
    enscount = 1
    if ensemble:
       enscount = algoinfo["enscount"]
    allsol = []
    modcurstates = {obstimes[timeindex]:deepcopy(curstates[timeindex]) for timeindex in xrange(len(obstimes))}
    node2slots = boundInferUtil.assignStateTimeSlots(modcurstates,smodel,G)
    node2vars = boundInferUtil.getNode2Vars(node2slots,smodel)
    linear,quad = boundInferUtil.estimateCoefs(node2vars,smodel,G,maxtime,APPROX,"normal")
    for index in xrange(enscount):
        if method == "greedy":
           history = greedyMaxGroupCoverage(deepcopy(node2vars),deepcopy(linear),deepcopy(quad),maxtime)
        elif method == "pipage":
           history = runLp(modcurstates,smodel,G,runfolder,deepcopy(node2slots),deepcopy(node2vars),deepcopy(linear),deepcopy(quad))
        history = HistoryUtil.lastFixSolution(history,G,smodel,"bound",infermode) 
        allsol.append(deepcopy(history)) 
    if ensemble:
       history = DisEns.getEnsembleSolution(allsol,infermode,curstates,G,smodel)
    if infermode == "Spreader":
       foundtime = None
       for time in sorted(allsol[0].keys()):
           if len(allsol[0][time][Trace.INFECTED]) != 0:
              foundtime = time
              break 
       history = allsol[0][foundtime], 0
       print "Inferred Spreader Solution at time {0} is: ".format(history[1],history[0])
    elif infermode == "History":     
       history = allsol[0]
       print "given trace info"
       for time in sorted(modcurstates.keys()):  
           print time,len(modcurstates[time][Trace.SUSCEPTIBLE]),len(modcurstates[time][Trace.EXPOSED]),len(modcurstates[time][Trace.INFECTED]),len(modcurstates[time][Trace.RECOVERED])        
       print "inferred history solution"
       for time in sorted(history.keys()):
           print time,len(history[time][Trace.SUSCEPTIBLE]),len(history[time][Trace.EXPOSED]),len(history[time][Trace.INFECTED]),len(history[time][Trace.RECOVERED])  
         
    if globals()["isTest"]:
       [QsapMinTest.checkVariableValidity(modcurstates,varname,smodel) for node in node2vars.keys() for varname in node2vars[node]]
       checkUtil.checkSolution(infermode,"bound",history,G,smodel,"notdetailed")   #min may not have to be at 0!
    InputOutput.writeHistoryResult(history,infermode,resultfile)   
        
            
if __name__ == "__main__":
    with gzip.open(sys.argv[1],"rb") as infile:  
        algopar = cPickle.load(infile)      
        curstates = cPickle.load(infile)
        obstimes = cPickle.load(infile)
        G = cPickle.load(infile)
        smodel = cPickle.load(infile)
        infermode = cPickle.load(infile)
        resultfile = cPickle.load(infile)
        runfolder = cPickle.load(infile)
    runner(algopar,curstates,obstimes,G,smodel,infermode,resultfile,runfolder)
