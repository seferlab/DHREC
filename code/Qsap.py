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
import numpy as np
import checkUtil
import cutUtil


class QsapTest():
     #Qsap Test Class
     def __init__():
         return

     @staticmethod
     def isSolutionSame(lpsol,flowsol):
         """compare 2 solutions and return true if same
         Args:
            lpsol:
            flowsol:
         Returns:   
            bool: true/false
         """
         if len(set(lpsol.keys()).symmetric_difference(flowsol.keys())) != 0:
            return False 
         for time in sorted(lpsol.keys()):
             for state in lpsol[time].keys():
                 if len(set(lpsol[time][state]).symmetric_difference(set(flowsol[time][state]))) != 0:
                    return False
         return True         
          
     @staticmethod
     def checkNodeVars(node2vars,smodel):
         """checks node var assignments
         Args:
            node2vars:
            smodel:
         Returns:
            bool: true/false   
         """
         for node in node2vars.keys():
             for varname in node2vars[node]:
                 times = HistoryUtil.var2Times(varname)
                 if smodel == "seir" and times.has_key(Trace.INFECTED) and not times.has_key(Trace.EXPOSED):
                    assert times[Trace.INFECTED] == 0
                 if smodel == "seir" and times.has_key(Trace.EXPOSED):
                    assert times[Trace.EXPOSED] >= 1    
                 if times.has_key(Trace.INFECTED) and times.has_key(Trace.EXPOSED):
                    assert times[Trace.INFECTED] > times[Trace.EXPOSED]
                 if times.has_key(Trace.RECOVERED) and times.has_key(Trace.EXPOSED):
                    assert times[Trace.RECOVERED] > times[Trace.EXPOSED]
                 if times.has_key(Trace.RECOVERED) and times.has_key(Trace.INFECTED):
                    assert times[Trace.RECOVERED] > times[Trace.INFECTED]
                 if smodel != "seir":
                    assert not times.has_key(Trace.EXPOSED)
                 if smodel == "si":
                    assert not times.has_key(Trace.RECOVERED) 
         return True
              
     @staticmethod
     def checkFlowResult(mappedG,node2terminal,node2labels,sortedlabels,node2assign,s,t):
         """checks flow result from boykov algo and assignment
         Args:
            mappedG:
            node2terminal:
            node2labels:
            sortedlabels:
            node2assign:
            s:
            t:
         Returns:
            bool : true/false
         """ 
         allnodes = set(node for node,label in node2terminal.keys())
         assert len(allnodes.symmetric_difference(node2labels.keys())) == 0
         snodes, tnodes = set(), set()
         for node in node2terminal.keys():
             if node2terminal[node] == s:
                snodes.add(node)
             elif node2terminal[node] == t:
                tnodes.add(node) 
         assert len(snodes) + len(tnodes) == mappedG.number_of_nodes() -2  and len(node2terminal.keys()) == mappedG.number_of_nodes() - 2
         for node in node2labels.keys():
             curlabel = sortedlabels[0]
             change = 0
             newnode = (node,curlabel)
             curterminal = node2terminal[newnode]
             for index in xrange(1,len(sortedlabels)):
                 label = sortedlabels[index]
                 newnode = (node,label)
                 if node2terminal[newnode] != curterminal:
                    curterminal = node2terminal[newnode]
                    change += 1
             assert change <= 1        
         for node in node2assign.keys():
             curlabel = node2assign[node]
             newnode = (node,curlabel)
             assert node2terminal[newnode] == s 
             curindex = sortedlabels.index(curlabel)
             for index in xrange(0,curindex+1):
                 label = sortedlabels[index]
                 newnode = (node,label)
                 assert node2terminal[newnode] == s
             for index in xrange(curindex+1,len(sortedlabels)):
                 label = sortedlabels[index]
                 newnode = (node,label)
                 assert node2terminal[newnode] == t    
         for node in node2assign.keys():
             assert node2assign[node] in node2labels[node]
             
         return True        
    
     @staticmethod
     def cutGraphCheck(cutG,sortedlabels,node2labels,G,quad,s,t,INFWEIGHT):
         """checks cut graph
         Args:  
            cutG:
            sortedlabels: 
            node2labels:
            G:
            quad:
            s:
            t:
            INFWEIGHT:
         Returns:
            bool: true/false
         """        
         assert cutG.number_of_nodes()-2 == len(sortedlabels) * len(node2labels.keys()) and cutG.number_of_selfloops() == 0          
         assert len(cutG.successors(t)) == 0 and len(cutG.predecessors(s)) == 0
         assert len(cutG.successors(s)) == len(node2labels.keys()) and len(cutG.predecessors(t)) == len(node2labels.keys())
         for node,label in cutG.successors(s):
             assert cutG[s][(node,label)]["weight"] == INFWEIGHT
         assert not cutG.has_edge(s,t) and not cutG.has_edge(t,s)
         curnodes = list(set(cutG.nodes()).difference(set([s,t])))
         for node1,label1 in curnodes:
             for node2,label2 in curnodes:
                 if node1 == node2:
                    continue
                 if label1 not in node2labels[node1] or label2 not in node2labels[node2]:
                    continue
                 if len(node2labels[node1]) == 1 and len(node2labels[node2]) == 1:
                    continue 
                 curval = 0.0
                 var1 = "x{0}?{1}".format(node1,label1)
                 var2 = "x{0}?{1}".format(node2,label2)
                 if quad.has_key((var2,var1)):
                    curval = quad[(var2,var1)]
                 elif quad.has_key((var1,var2)):
                    curval = quad[(var1,var2)]   
                 in1 = sortedlabels.index(label1)
                 in2 = sortedlabels.index(label2)
                 cursum = 0.0 
                 for first,second,sendnode,recnode in [(in1+1,in2+1,node1,node2),(in2+1,in1+1,node2,node1)]:
                     for index1 in xrange(0,first):
                         sendlabel = sortedlabels[index1]
                         newsendnode = (sendnode,sendlabel)
                         for index2 in xrange(second,len(sortedlabels)):
                             reclabel = sortedlabels[index2]
                             newrecnode = (recnode,reclabel)
                             if cutG.has_edge(newsendnode,newrecnode):
                                cursum += cutG[newsendnode][newrecnode]["weight"]
                 assert abs(curval - cursum) <= 0.0001     
                              
         for gennode1,gennode2 in cutG.edges():
             assert cutG[gennode1][gennode2]["weight"] >= 0
             if gennode1 == s:
                node2,label2 = gennode2
                assert cutG[gennode1][gennode2]["weight"] == INFWEIGHT
             elif gennode2 == t:
                node1,label1 = gennode1
                if label1 in node2labels[node1]:
                   assert cutG[gennode1][gennode2]["weight"] != INFWEIGHT
                else:   
                   assert cutG[gennode1][gennode2]["weight"] == INFWEIGHT
             else:
                node1,label1 = gennode1
                node2,label2 = gennode2
                if node1 != node2:
                   assert (G.has_edge(node1,node2) or G.has_edge(node2,node1)) and cutG[gennode1][gennode2] != INFWEIGHT
                else: 
                   preindex = sortedlabels.index(label2)
                   curindex = sortedlabels.index(label1) 
                   assert abs(curindex-preindex) == 1 
                   if (curindex - preindex) == 1:
                      assert cutG[gennode1][gennode2]["weight"] == INFWEIGHT
                   elif (curindex - preindex) == -1:
                      if label1 in node2labels[node1]:
                         assert cutG[gennode1][gennode2]["weight"] != INFWEIGHT
                      else:
                         assert cutG[gennode1][gennode2]["weight"] == INFWEIGHT
         return True
           
     @staticmethod
     def transDependencyCheck(labelG,status):
        """transitive check in label graph
        Args:
           labelG:
           status:
        Returns:
           bool: true/false 
        """
        ties = []
        for in1 in xrange(len(labelG.nodes())):
            node1 = labelG.nodes()[in1]
            for in2 in xrange(in1+1,len(labelG.nodes())):
                node2 = labelG.nodes()[in2]
                if not labelG.has_edge(node1,node2) and not labelG.has_edge(node2,node1):
                   ties.append((node1,node2)) 
        print ties
        if status == "weak":
           return True
        elif status == "normal":
           intrans = 0
           for node1,node2,node3 in itertools.product(labelG.nodes(),labelG.nodes(),labelG.nodes()):
               if labelG.has_edge(node1,node2) and labelG.has_edge(node2,node3) and not labelG.has_edge(node1,node3):
                  intrans += 1
           return intrans == 0
     
     @staticmethod
     def coefCheck(quad,linear,smodel,G):
         """coef check 
         Args:
           quad:
           linear: 
           smodel:
           G:
         Returns:
           bool: true/false
         """
         if smodel in ["si","sir"]:
            sstate,rstate = Trace.INFECTED, Trace.INFECTED
         elif smodel == "seir":
            sstate,rstate =  Trace.INFECTED, Trace.EXPOSED    
         curquadkeys = quad.keys() 
         for var1,var2 in curquadkeys:
             times1 = HistoryUtil.var2Times(var1)
             times2 = HistoryUtil.var2Times(var2)
             node1 = int(var1.split("?")[0].replace("x",""))
             node2 = int(var2.split("?")[0].replace("x",""))
             assert HistoryUtil.checkAffect(times2,times1,smodel,sstate,rstate,"weak") and G.has_edge(node2,node1) 
             assert not quad.has_key((var2,var1)) 
             assert quad[(var1,var2)] > 0    
         for varname in linear.keys():
             assert linear[varname] > 0    
         if smodel in ["si","sir"]:
            for var1,var2 in quad.keys():
                times1 = HistoryUtil.var2Times(var1)
                times2 = HistoryUtil.var2Times(var2)
                assert times1[Trace.INFECTED] > times2[Trace.INFECTED] 
         return True        
     
     @staticmethod
     def checkGraphValidness(G,sol,smodel):
         """checks graph validness
         Args:
           G:
           sol: sol in times format
           smodel:
         Returns:
           bool:
         """
         #if smodel in ["sir","seir"]:
         recovtimes = {node:time  for time in sol.keys() for node in sol[time][Trace.RECOVERED]}
         rectimes,sendtimes = {},{}
         if smodel in ["si","sir"]:
            rectimes = {time:set(sol[time][Trace.INFECTED]) for time in sol.keys()}
            sendtimes = deepcopy(rectimes)
         elif smodel == "seir":
            rectimes = {time:set(sol[time][Trace.EXPOSED]) for time in sol.keys()}
            sendtimes = {time:set(sol[time][Trace.INFECTED]) for time in sol.keys()}
         sortedtimes = sorted(sol.keys())
         for timeindex in xrange(1,len(sortedtimes)):
             time = sortedtimes[timeindex]
             curnodes = set(rectimes[time])
             for timeindex2 in xrange(0,timeindex):
                 time2 = sortedtimes[timeindex2]
                 curnodes2 = set(sendtimes[time2])
                 delnodes = set()
                 for node in curnodes:
                     flag = False
                     for node2 in curnodes2:
                         if G.has_edge(node2,node) and (not recovtimes.has_key(node2) or (recovtimes.has_key(node2) and recovtimes[node2] >= time)):
                            flag = True    
                            break
                     if flag:
                        delnodes.add(node)
                 curnodes = curnodes.difference(delnodes)
             #assert len(curnodes) == 0
             print time,len(curnodes),len(rectimes[time])  
         return True    
         
     @staticmethod
     def checkVariableValidity(modcurstates,varname,smodel):
         """check whether variable is valid on given modcurstates
         Args:
           modcurstates:
           varname:
           smodel:
         """
         node = int(varname.replace("x","").split("?")[0])
         sortedtimes = sorted(modcurstates.keys())
         nodetimes = HistoryUtil.var2Times(varname)
         if smodel == "si":
            for time in sortedtimes:
                if time < nodetimes[Trace.INFECTED]:
                   assert node in modcurstates[time][Trace.SUSCEPTIBLE]
                else:
                   assert node in modcurstates[time][Trace.INFECTED] 
         elif smodel == "sir":
            if nodetimes.has_key(Trace.INFECTED) and nodetimes.has_key(Trace.RECOVERED):
               for time in sortedtimes:
                   if time < nodetimes[Trace.INFECTED]:
                      assert node in modcurstates[time][Trace.SUSCEPTIBLE]
                   elif time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.INFECTED]
                   elif time >= nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.RECOVERED]    
            elif nodetimes.has_key(Trace.INFECTED):
               for time in sortedtimes:
                   if time < nodetimes[Trace.INFECTED]:
                      assert node in modcurstates[time][Trace.SUSCEPTIBLE]
                   else:
                      assert node in modcurstates[time][Trace.INFECTED]        
         elif smodel == "seir":
            if nodetimes.has_key(Trace.EXPOSED) and nodetimes.has_key(Trace.INFECTED) and nodetimes.has_key(Trace.RECOVERED):
               assert nodetimes[Trace.EXPOSED] > 0 
               for time in sortedtimes:
                   if time < nodetimes[Trace.EXPOSED]:
                      assert node in modcurstates[time][Trace.SUSCEPTIBLE]
                   elif time >= nodetimes[Trace.EXPOSED] and time < nodetimes[Trace.INFECTED]:
                      assert node in modcurstates[time][Trace.EXPOSED]
                   elif time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.INFECTED]   
                   elif time >= nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.RECOVERED]   
            elif nodetimes.has_key(Trace.INFECTED) and nodetimes.has_key(Trace.RECOVERED):
               assert nodetimes[Trace.INFECTED] == 0  
               for time in sortedtimes:
                   if time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.INFECTED]
                   elif time >= nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.RECOVERED]       
            elif nodetimes.has_key(Trace.EXPOSED) and nodetimes.has_key(Trace.INFECTED):
               assert nodetimes[Trace.EXPOSED] > 0 
               for time in sortedtimes:
                   if time < nodetimes[Trace.EXPOSED]:
                      assert node in modcurstates[time][Trace.SUSCEPTIBLE]
                   elif time >= nodetimes[Trace.EXPOSED] and time < nodetimes[Trace.INFECTED]:
                      assert node in modcurstates[time][Trace.EXPOSED]
                   elif time >= nodetimes[Trace.INFECTED]:
                      assert node in modcurstates[time][Trace.INFECTED] 
            elif nodetimes.has_key(Trace.INFECTED):
               assert nodetimes[Trace.INFECTED] == 0 
            elif nodetimes.has_key(Trace.EXPOSED):
               assert nodetimes[Trace.EXPOSED] > 0 
         return True
      
     @staticmethod 
     def lpSolutionCheck(paraminfo,solinfo):
         """checks validness of cplex solution
         Args: 
            paraminfo:
            solinfo: 
         Returns:
            bool: true/false
         """
         [node2slots,node2vars,modcurstates,quad,linear,G,smodel] = paraminfo
         [retvalues,sol,lpobjval,retobjval,isfrac,node2label] = solinfo
         assert node2label != None
         states = [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
         sortedtimes = sorted(modcurstates.keys())
         seennodes = set([node for state in states for node in modcurstates[max(sortedtimes)][state]])
         snaphist = HistoryUtil.buildSnapshotFromTransTimes(sol,smodel,G)
         for time in snaphist.keys():
             assert sum([len(snaphist[time][state]) for state in snaphist[time].keys()]) == G.number_of_nodes()
         for time in sortedtimes:
             posiset = set(modcurstates[time][Trace.INFECTED]).union(set(modcurstates[time][Trace.RECOVERED]))
             for snaptime in sorted(snaphist.keys()):
                 if snaptime > time:
                    break
                 for node in snaphist[snaptime][Trace.INFECTED]:
                     if node in posiset:
                        posiset.remove(node)
             assert len(posiset) == 0                         
             
         assert len(seennodes.difference(set(node2slots.keys()))) == 0 and len(set(node2slots.keys()).difference(seennodes)) == 0  
         assert len(set(node2slots.keys()).intersection(set(node2vars.keys()))) == len(set(node2slots.keys())) and len(set(node2slots.keys()).intersection(seennodes)) == len(seennodes)  
         for node in node2vars.keys():
             assert len(node2vars[node]) > 0
         singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
         doublevars = [var for var in retvalues.keys() if var.find("_") != -1]
         assert len(singlevars) != 0
         for doublevar in doublevars:
             var1,var2 = doublevar.split("_")
             newvar = "{0}_{1}".format(var2,var1)
             assert newvar in doublevars
             assert retvalues[doublevar] == retvalues[newvar]             
         for doublevar in doublevars:
             var1,var2 = doublevar.split("_")
             assert var1 in singlevars or var2 in singlevars
         startnodes = set()
         for node in node2label.keys():    
             varname = "x{0}?{1}".format(node,node2label[node])
             temptimes = HistoryUtil.var2Times(varname)
             if temptimes.has_key(Trace.INFECTED) and temptimes[Trace.INFECTED] == 0:
                startnodes.add(node)                 
         difset = startnodes.symmetric_difference(set(sol[0][Trace.INFECTED]))
         for item in difset:
             varname = "x{0}?{1}".format(item,node2label[item])
             temptimes = HistoryUtil.var2Times(varname)
             print item,temptimes,varname      
         assert len(startnodes.symmetric_difference(set(sol[0][Trace.INFECTED]))) == 0
         estobjval = reestimateObjVal(retvalues,linear,quad)
         assert abs(estobjval - lpobjval) <= 0.1  
         
         if not isfrac:
            for varname in singlevars:
                assert retvalues[varname] == 1.0
            for varname in doublevars:
                assert retvalues[varname] == 1.0  
            for node in node2vars.keys():
                seenset = set()
                for myslim in node2vars[node]:    
                    if retvalues.has_key(myslim): 
                       seenset.add(myslim)
                assert len(seenset) == 1    
            node2val = {}
            for varname in singlevars:
                node = int(varname.replace("x","").split("?")[0])
                node2val.setdefault(node,{})
                node2val[node][varname] = retvalues[varname] 
            for node in node2val.keys():
                assert len(node2val[node].keys()) == 1 and node2val[node].values()[0] == 1.0 
                assert sum([node2val[node][varname] for varname in node2val[node].keys()]) == 1.0  
            allnodes = [int(singlevar.replace("x","").split("?")[0]) for singlevar in singlevars]
            assert len(allnodes) == len(set(allnodes))
            
         againobjval = estimateObjVal(node2label,linear,quad,G)
         assert abs(againobjval - retobjval) <= 0.1  
         
         return True

     
def reestimateObjVal(retvalues,linear,quad):
    """estimate obj val from retvalues
    Args:
       retvalues:
       linear:
       quad:
    Returns:
       estobjval:   
    """         
    singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
    doublevars = [var for var in retvalues.keys() if var.find("_") != -1]
    estobjval = sum([retvalues[singlevar] * linear[singlevar] for singlevar in singlevars if linear.has_key(singlevar)])
    for doublevar in doublevars:
        var1,var2 = doublevar.split("_")
        if quad.has_key((var1,var2)):
           estobjval += quad[(var1,var2)] * retvalues[doublevar]  
    return estobjval     

     
def getTime(sortedtimes,searchval):
    """returns the previous time(or 0 if searchval is at the beginning of sortedtimes)
    Args:
       sortedtimes:
       searchval:
    Returns:
       time:
    """
    tindex = sortedtimes.index(searchval) 
    if tindex != 0:
       return sortedtimes[tindex-1] + 1
    return 0

def getLastTime(sortedtimes,modcurstates,node2slots,node,state):                  
    """gets Last possible time
    Args:
       sortedtimes:
       modcurstates:
       node2slots:
       node:
       state:
    Returns:
       time:
    """
    if state == Trace.SUSCEPTIBLE:
       inittime = -1
       for time in sortedtimes:
           if node in modcurstates[time][Trace.SUSCEPTIBLE]:
              inittime = time
           else:     
              break
       return inittime     
    
    prein = sortedtimes.index(node2slots[node][state]["end"]) 
    time = node2slots[node][state]["end"]
    for timein in xrange(prein+1,len(sortedtimes)):
        if node in modcurstates[sortedtimes[timein]][state]:
           time = sortedtimes[timein]
        else:
           break 
    return time

def assignStateTimeSlots(modcurstates,smodel,G):
    """assign time slots for nodes and states
    Args:
       modcurstates:
       smodel:
       G:
    Returns:
       node2slots: 
    """
    allnodes = set(node for state in modcurstates.values()[0].keys() for node in modcurstates.values()[0][state]) 
    sortedtimes = sorted(modcurstates.keys()) 
    states = [Trace.INFECTED,Trace.EXPOSED,Trace.RECOVERED]
    node2slots = {node: {state: {"end":None} for state in states} for node in allnodes}
    for time in sortedtimes:
        for state in modcurstates[time].keys():
            if state == Trace.SUSCEPTIBLE:
               continue
            for node in modcurstates[time][state]:
                if node2slots[node][state]["end"] == None:
                   node2slots[node][state]["end"] = time 
    for node in node2slots.keys():
        for state in node2slots[node].keys():
            if node2slots[node][state]["end"] == None:
               del node2slots[node][state]  
        if len(node2slots[node].keys()) == 0:
           del node2slots[node]

    if globals()["isTest"]:     
       seennodes = set([node for state in states for node in modcurstates[max(sortedtimes)][state]])
       assert len(seennodes.difference(set(node2slots.keys()))) == 0 and len(set(node2slots.keys()).difference(seennodes)) == 0    
       
    for node in node2slots.keys(): 
        if smodel == "si":
           node2slots[node][Trace.INFECTED]["begin"] = getTime(sortedtimes,node2slots[node][Trace.INFECTED]["end"])
        elif smodel == "sir":    
           if not node2slots[node].has_key(Trace.INFECTED) and node2slots[node].has_key(Trace.RECOVERED):
              pretime = getTime(sortedtimes,node2slots[node][Trace.RECOVERED]["end"])
              node2slots[node] = {Trace.INFECTED: {"begin": pretime, "end": node2slots[node][Trace.RECOVERED]["end"]-1}, Trace.RECOVERED: {"begin": pretime+1, "end": node2slots[node][Trace.RECOVERED]["end"]}}
           elif node2slots[node].has_key(Trace.INFECTED):
              node2slots[node][Trace.INFECTED]["begin"] = getTime(sortedtimes,node2slots[node][Trace.INFECTED]["end"])           
              if node2slots[node].has_key(Trace.RECOVERED): 
                 node2slots[node][Trace.RECOVERED]["begin"] = getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.INFECTED) + 1
        elif smodel == "seir":
           if node2slots[node].has_key(Trace.EXPOSED) and node2slots[node].has_key(Trace.RECOVERED) and node2slots[node].has_key(Trace.INFECTED):
              node2slots[node][Trace.EXPOSED]["begin"] = max(1,getTime(sortedtimes,node2slots[node][Trace.EXPOSED]["end"]))
              node2slots[node][Trace.INFECTED]["begin"] = getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.EXPOSED) + 1
              node2slots[node][Trace.RECOVERED]["begin"] = getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.INFECTED) + 1
           elif node2slots[node].has_key(Trace.EXPOSED) and node2slots[node].has_key(Trace.RECOVERED):
              node2slots[node][Trace.EXPOSED]["begin"] = max(1,getTime(sortedtimes,node2slots[node][Trace.EXPOSED]["end"]))
              node2slots[node][Trace.INFECTED] = {}
              node2slots[node][Trace.INFECTED]["begin"] = getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.EXPOSED) + 1
              node2slots[node][Trace.INFECTED]["end"] = node2slots[node][Trace.RECOVERED]["end"] - 1
              node2slots[node][Trace.RECOVERED]["begin"] = getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.EXPOSED) + 1
           elif node2slots[node].has_key(Trace.EXPOSED) and node2slots[node].has_key(Trace.INFECTED):
              node2slots[node][Trace.EXPOSED]["begin"] = max(1,getTime(sortedtimes,node2slots[node][Trace.EXPOSED]["end"]))
              node2slots[node][Trace.INFECTED]["begin"] = getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.EXPOSED) + 1
           elif node2slots[node].has_key(Trace.INFECTED) and node2slots[node].has_key(Trace.RECOVERED): #?
              node2slots[node][Trace.EXPOSED]= {"begin":max(1,getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.SUSCEPTIBLE)+1), "end":node2slots[node][Trace.INFECTED]["end"] - 1}
              node2slots[node][Trace.INFECTED]["begin"] = getTime(sortedtimes,node2slots[node][Trace.INFECTED]["end"])
              node2slots[node][Trace.RECOVERED]["begin"] = getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.INFECTED) + 1
           elif node2slots[node].has_key(Trace.RECOVERED): #?
              node2slots[node][Trace.RECOVERED]["begin"] = getTime(sortedtimes,node2slots[node][Trace.RECOVERED]["end"])
              node2slots[node][Trace.EXPOSED] = {"begin":max(1,getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.SUSCEPTIBLE)+1),"end":node2slots[node][Trace.RECOVERED]["end"] - 2}
              node2slots[node][Trace.INFECTED] = {"begin":max(0,getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.SUSCEPTIBLE)+1),"end":node2slots[node][Trace.RECOVERED]["end"] - 1}
           elif node2slots[node].has_key(Trace.INFECTED): #?
              node2slots[node][Trace.INFECTED]["begin"] = getTime(sortedtimes,node2slots[node][Trace.INFECTED]["end"])
              node2slots[node][Trace.EXPOSED] = {"begin":max(1,getLastTime(sortedtimes,modcurstates,node2slots,node,Trace.SUSCEPTIBLE)+1),"end":node2slots[node][Trace.INFECTED]["end"] - 1}  
           elif node2slots[node].has_key(Trace.EXPOSED):
              node2slots[node][Trace.EXPOSED]["begin"] = max(1,getTime(sortedtimes,node2slots[node][Trace.EXPOSED]["end"])) 
        if len(G.predecessors(node)) == 0:
           node2slots[node][Trace.INFECTED] = {"begin":0,"end":0}
    return node2slots 


def getNode2Vars(node2slots,smodel):
    """gets node2vars
    Args:
       node2slots:
       smodel:
    Returns:
       node2vars:
    """
    node2vars = {}
    for node in node2slots.keys():
        if smodel == "si":
           for time in xrange(node2slots[node][Trace.INFECTED]["begin"],node2slots[node][Trace.INFECTED]["end"]+1):
               varname = "x{0}?{1}?{2}".format(node,Trace.INFECTED,time)
               node2vars.setdefault(node,set()).add(varname)
        elif smodel == "sir":
           if node2slots[node].has_key(Trace.INFECTED) and node2slots[node].has_key(Trace.RECOVERED):
              for itime in xrange(node2slots[node][Trace.INFECTED]["begin"],node2slots[node][Trace.INFECTED]["end"]+1):
                  for rtime in xrange(node2slots[node][Trace.RECOVERED]["begin"],node2slots[node][Trace.RECOVERED]["end"]+1):
                      if itime < rtime:
                         varname = "x{0}?{1}?{2}?{3}?{4}".format(node,Trace.INFECTED,itime,Trace.RECOVERED,rtime)
                         node2vars.setdefault(node,set()).add(varname)
           elif not node2slots[node].has_key(Trace.RECOVERED):   
              assert node2slots[node].has_key(Trace.INFECTED)
              for itime in xrange(node2slots[node][Trace.INFECTED]["begin"],node2slots[node][Trace.INFECTED]["end"]+1):
                  varname = "x{0}?{1}?{2}".format(node,Trace.INFECTED,itime)
                  node2vars.setdefault(node,set()).add(varname)     
        elif smodel == "seir": 
           if node2slots[node].has_key(Trace.EXPOSED) and node2slots[node].has_key(Trace.INFECTED) and node2slots[node].has_key(Trace.RECOVERED):
              for etime in xrange(node2slots[node][Trace.EXPOSED]["begin"],node2slots[node][Trace.EXPOSED]["end"]+1):
                  for itime in xrange(node2slots[node][Trace.INFECTED]["begin"],node2slots[node][Trace.INFECTED]["end"]+1):
                      if etime >= itime:
                         continue 
                      for rtime in xrange(node2slots[node][Trace.RECOVERED]["begin"],node2slots[node][Trace.RECOVERED]["end"]+1):
                          if itime < rtime:
                             varname = "x{0}?{1}?{2}?{3}?{4}?{5}?{6}".format(node,Trace.EXPOSED,etime,Trace.INFECTED,itime,Trace.RECOVERED,rtime)
                             node2vars.setdefault(node,set()).add(varname)
              if node2slots[node][Trace.INFECTED]["begin"] == 0:
                 for rtime in xrange(node2slots[node][Trace.RECOVERED]["begin"],node2slots[node][Trace.RECOVERED]["end"]+1):
                     if rtime > 0 :
                        varname = "x{0}?{1}?{2}?{3}?{4}".format(node,Trace.INFECTED,0,Trace.RECOVERED,rtime)
                        node2vars.setdefault(node,set()).add(varname)
           elif node2slots[node].has_key(Trace.EXPOSED) and node2slots[node].has_key(Trace.INFECTED): 
              for etime in xrange(node2slots[node][Trace.EXPOSED]["begin"],node2slots[node][Trace.EXPOSED]["end"]+1):
                  for itime in xrange(node2slots[node][Trace.INFECTED]["begin"],node2slots[node][Trace.INFECTED]["end"]+1):
                      if etime >= itime:
                         continue 
                      varname = "x{0}?{1}?{2}?{3}?{4}".format(node,Trace.EXPOSED,etime,Trace.INFECTED,itime)
                      node2vars.setdefault(node,set()).add(varname)
              if node2slots[node][Trace.INFECTED]["begin"] == 0:
                 varname = "x{0}?{1}?{2}".format(node,Trace.INFECTED,0)
                 node2vars.setdefault(node,set()).add(varname)
           elif node2slots[node].has_key(Trace.EXPOSED):
              for etime in xrange(node2slots[node][Trace.EXPOSED]["begin"],node2slots[node][Trace.EXPOSED]["end"]+1):
                  varname = "x{0}?{1}?{2}".format(node,Trace.EXPOSED,etime)
                  node2vars.setdefault(node,set()).add(varname)
                  
    if globals()["isTest"]:
       QsapTest.checkNodeVars(node2vars,smodel)
    return node2vars

         
def greedyMaxGroupCoverage(node2vars,linear,quad,nodevarG,maxtime):
    """implements greedy max group coverage algorithm for max case (1/2 approximation)
    Args:
       node2vars:
       linear:
       quad:
       nodevarG:
       maxtime:
    Returns:
       history:
    """
    cursol = set((node,node2vars[node][0]) for node in node2vars.keys() if len(node2vars[node]) == 1)
    curnodes =[node for node in node2vars.keys() if len(node2vars[node]) > 1]
    optval = sum([linear[varname] for node,varname in cursol])
    optval += sum([quad[varname][prenodevar] for node,varname in cursol for prenodevar in set(nodevarG.predecessors((node,varname)))])
    optval -= sum([quad[nodevar1][nodevar2] for nodevar1 in cursol for nodevar2 in cursol if quad.has_key((nodevar1,nodevar2))])
    while len(curnodes) != 0:
        maxnodevar,maxinc = None, 0.0
        for node in curnodes:
            for varname in node2vars[node]:
                inc = linear[varname]
                inc += sum([quad[varname][prenodevar] for prenodevar in set(nodevarG.predecessors((node,varname)))])
                inc -= sum([quad[nodevar1][nodevar2] for nodevar1 in cursol for nodevar2 in cursol if quad.has_key((nodevar1,nodevar2))])
                assert inc >= 0
                if inc >= maxinc:
                   maxnodevar = (node,varname)
                   maxinc = inc  
        assert maxinc >= 0        
        assert maxnodevar != None
        curnodes.remove(maxnodevar[0])
        cursol.add(maxnodevar)
        optval += inc    
    print cursol
    print optval
    history = {time: {state: set() for state in [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]} for time in xrange(maxtime+1)}
    print history
    exit(1)
    return history




def estimateObjVal(node2assign,linear,quad,G):
    """estimate objective value
    Args:
      node2assign:
      linear:
      quad:
      G:
    Returns:
      objval = 0.0 
    """    
    objval = 0.0   
    for node in node2assign.keys():
        label = node2assign[node]
        varname = "x{0}?{1}".format(node,label)
        if linear.has_key(varname):
           objval += linear[varname]
    curnodes = list(node2assign.keys())  
    for node1 in curnodes:
        label1 = node2assign[node1]
        var1 = "x{0}?{1}".format(node1,label1)
        for node2 in curnodes:
            label2 = node2assign[node2]
            var2 = "x{0}?{1}".format(node2,label2)
            if G.has_edge(node2,node1) and quad.has_key((var1,var2)):
               objval += quad[(var1,var2)]    
    return objval
           

def estimateTransCoefs(node,times1,varname1,smodel,sstate,rstate,transdist,G,node2vars,ZERO,APPROX):
    """estimate transition coefs(s2i or s2e)
    Args:
       node:
       times1:
       varname1:
       smodel:
       sstate:
       rstate:
       transdist:
       G:
       node2vars:
       ZERO:
       APPROX:
    Returns:
       quad
    """
    quad = {} 
    for neighnode in set(G.predecessors(node)).intersection(set(node2vars.keys())):
        for varname2 in node2vars[neighnode]:
            times2 = HistoryUtil.var2Times(varname2)
            if HistoryUtil.checkAffect(times2,times1,smodel,sstate,rstate,"weak"):
               f1 = Dist.genPartDist(G[neighnode][node][transdist][0],G[neighnode][node][transdist][1],"reverseCdf",G[neighnode][node][Trace.SPROB])
               fdist = Dist.genPdf(f1)
               fdist[0] = 1.0
               p1 = Dist.genPartDist(G[neighnode][node][transdist][0],G[neighnode][node][transdist][1],"normal",G[neighnode][node][Trace.SPROB])
               pdist = Dist.genPdf(p1)
               diftime = times1[rstate] - times2[sstate]
               assert diftime > 0
               if APPROX == "reliability": 
                  val = ZERO
                  if fdist.has_key(diftime-1) and fdist[diftime-1] > ZERO:
                     val = fdist[diftime-1] 
               elif APPROX == "reliability2": 
                  val = ZERO
                  if fdist.has_key(diftime) and fdist[diftime] > ZERO:
                     val = fdist[diftime]        
               elif APPROX == "prob": 
                  val = ZERO
                  if pdist.has_key(diftime) and pdist[diftime] > ZERO:
                     val = pdist[diftime]
               elif APPROX == "ratio":
                  part1,part2 = ZERO,ZERO
                  if fdist.has_key(diftime) and fdist[diftime] > ZERO:
                     part1 = fdist[diftime]
                  if fdist.has_key(diftime-1) and fdist[diftime-1] > ZERO:
                     part2 = fdist[diftime-1]
                  val = part2 / part1
                  assert val >= 1.0
                  val = 1.0 / val 
               if val != 1.0:
                  quad[(varname1,varname2)] = -1.0*math.log(val) 
    return quad


def estimateEndoCoefs(node,times,G,transtype,ZERO,maxtime):
    """estimate endogenous transition coefs(i2r,e2i) 
    Args:
       node:
       times:
       G:
       transtype:
       ZERO:
       maxtime:
    Returns:
       coef:
    """       
    assert transtype in ["i2r","e2i"]
    statemap = {"s":Trace.SUSCEPTIBLE,"e":Trace.EXPOSED,"i":Trace.INFECTED,"r":Trace.RECOVERED}
    curstate,nextstate = [statemap[state] for state in transtype.split("2")]
    if transtype == "i2r":
       transname = Trace.I2R  
    elif transtype == "e2i":
       transname = Trace.E2I  

    coef = 0.0
    if times.has_key(curstate) and times.has_key(nextstate):
       p1 = Dist.genPartDist(G.node[node][transname][0],G.node[node][transname][1],"normal")
       pdist = Dist.genPdf(p1)  
       diftime = times[nextstate] - times[curstate]
       if pdist.has_key(diftime) and pdist[diftime] >= ZERO:
          coef = -1.0 * math.log(pdist[diftime])
       else:
          coef = -1.0 * math.log(ZERO)   
    elif times.has_key(curstate):
       p1 = Dist.genPartDist(G.node[node][transname][0],G.node[node][transname][1],"reverseCdf")
       fdist = Dist.genPdf(p1) 
       diftime = maxtime - times[curstate]
       if diftime == 0:
          coef = 0.0
       elif fdist.has_key(diftime) and fdist[diftime] >= ZERO:
          coef = -1.0 * math.log(fdist[diftime])
       else:
          coef = -1.0 * math.log(ZERO)
    return coef

                     
def estimateNoneTransCoef(snodes,G,node2vars,smodel,ZERO,maxtime):
    """estimates none trans coefs
    Args:
       snodes:
       G:
       node2vars:
       smodel:
       ZERO:
       maxtime:
    Returns:
       linear:
    """
    if smodel in ["si","sir"]:
       transdist = Trace.S2I
    elif smodel == "seir":
       transdist = Trace.S2E   
    linear = {}
    for node in snodes:
        for neighnode in set(G.predecessors(node)).intersection(set(node2vars.keys())):
            for varname in node2vars[neighnode]:
                mytimes = HistoryUtil.var2Times(varname)
                if not mytimes.has_key(Trace.INFECTED):
                   continue 
                p1 = Dist.genPartDist(G[neighnode][node][transdist][0],G[neighnode][node][transdist][1],"reverseCdf")
                fdist = Dist.genPdf(p1)
                if not mytimes.has_key(Trace.RECOVERED):
                   diftime = maxtime - mytimes[Trace.INFECTED] 
                else:
                   diftime = mytimes[Trace.RECOVERED] - mytimes[Trace.INFECTED]
                if diftime == 0:
                   continue    
                linear.setdefault(varname,0.0) 
                if fdist.has_key(diftime) and fdist[diftime] >= ZERO:
                   linear[varname] += -1.0 * math.log(fdist[diftime])
                else:
                   linear[varname] += -1.0 * math.log(ZERO)
    return linear
                     

def estimateCoefs(node2vars,smodel,G,maxtime,APPROX):
    """estimate coefficients
    Args:
      node2vars:
      smodel:
      G:
      maxtime:
      APPROX:
    Returns:
      linear:
      quad: 
    """
    ZERO = 0.000000000000000001 
    if smodel in ["si","sir"]:
       sstate,rstate = Trace.INFECTED,Trace.INFECTED
       transdist = Trace.S2I
    elif smodel == "seir":
       sstate,rstate = Trace.INFECTED,Trace.EXPOSED
       transdist = Trace.S2E   
    linear = {varname:0.0 for node in node2vars.keys() for varname in node2vars[node]}
    quad = {} 
    snodes = set(G.nodes()).difference(set(node2vars.keys()))
    for node in node2vars.keys():
        for varname1 in node2vars[node]:
            times1 = HistoryUtil.var2Times(varname1)
            localquad = estimateTransCoefs(node,times1,varname1,smodel,sstate,rstate,transdist,G,node2vars,ZERO,APPROX)
            for curvarname1,varname2 in localquad.keys(): 
                assert not quad.has_key((varname1,varname2)) 
                quad[(varname1,varname2)] = localquad[(varname1,varname2)]
            if smodel in ["sir","seir"]:
               coef = estimateEndoCoefs(node,times1,G,"i2r",ZERO,maxtime)
               if coef != 0.0:
                  linear[varname1] += coef 
            if smodel == "seir":
               coef = estimateEndoCoefs(node,times1,G,"e2i",ZERO,maxtime)
               if coef != 0.0:
                  linear[varname1] += coef             
    locallinear = estimateNoneTransCoef(snodes,G,node2vars,smodel,ZERO,maxtime)
    for key in locallinear.keys():
        linear[key] += locallinear[key]
    for key in linear.keys():
        if linear[key] == 0.0 or linear[key] == -0.0:
           del linear[key] 
    if APPROX == "ratio":
       applinear = {}
       for node in node2vars.keys():
           for mainvar in node2vars[node]:
               if not linear.has_key(mainvar):
                  continue 
               for sidevar in set(node2vars[node]).difference(set(mainvar)):
                   applinear.setdefault(sidevar,0.0)
                   applinear[sidevar] += linear[mainvar]
       linear = deepcopy(applinear)  
    return linear,quad
 

def genMLBoundConst(linear,quad,node2vars,G):
    """generates Metric Labeling bound constraints
    Args:
       linear:
       quad:
       node2vars:
       G:
    Returns:
       boundstr:   
    """
    boundstr = "Bounds\n"
    boundstr += "\n".join(["0 <= {0} <= 1.0".format(varname) for varname in linear.keys()]) + "\n"
    #for var1,var2 in quad.keys():
    #    varname = "{0}_{1}".format(var1,var2)
    #    boundstr += "0 <= {0} <= 1.0\n".format(varname)
    for node in node2vars.keys():
        for var1 in node2vars[node]:
            tnodes = set(G.predecessors(node)).union(set(G.successors(node)))
            affectnodes = tnodes.intersection(node2vars.keys())
            for neighnode in affectnodes:
                for var2 in node2vars[neighnode]:
                    boundstr += "0 <= {0}_{1} <= 1.0\n".format(var1,var2)
                    boundstr += "0 <= {0}_{1} <= 1.0\n".format(var2,var1)
    return boundstr


def genMLObj(linear,quad,algoinfo):
    """generates QSAP objective function 
    Args:
      linear:
      quad:
      algoinfo:
    Returns:
      objstr:
    """   
    APPROX = algoinfo["approx"]
    isPlus = lambda x: "+" if x >= 0 else " "
    if APPROX == "ratio":
       objstr = "Maximize\n obj: "
    else:    
       objstr = "Minimize\n obj: "
    objstr += " ".join([" {0} %.10f {1} ".format(isPlus(linear[var1]),var1) %linear[var1]  for var1 in linear.keys()])       
    objstr += " ".join([" {0} %.10f {1} ".format(isPlus(quad[(var1,var2)]),var1+"_"+var2) %quad[(var1,var2)] for var1,var2 in quad.keys()])
    return objstr
       
def genMLMainConst(node2vars,linear,smodel,G,quad):
    """generates constraints 
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
    
    for node in node2vars.keys():
        for varname in node2vars[node]:
            tnodes = set(G.predecessors(node)).union(set(G.successors(node)))
            affectnodes = tnodes.intersection(node2vars.keys())
            #affectnodes = set(G.predecessors(node)).intersection(node2vars.keys())
            if len(affectnodes) != 0:
               for neighnode in affectnodes:
                   putvars = [" 1.0 {0}_{1} ".format(varname,var2) for var2 in node2vars[neighnode]]
                   consstr += " + ".join(putvars) + " - 1.0 {0} = 0.0\n".format(varname)
                   for var2 in node2vars[neighnode]:
                       consstr += " {0}_{1} - {1}_{0} = 0\n".format(var2,varname)
    
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
    

def runCode(consstr,objstr,boundstr,runfolder):
    """runs metric labeling code
    Args:
        consstr: constraint string
        objstr: objective function string
        boundstr: boundary string
        runfolder: run folder
    Returns:
        edge2val : edge value mapping
    """
    PNUM = 1
    algoparstr = "qsap"
    outlppath = "{0}/{1}.lp".format(runfolder,algoparstr)
    with open(outlppath,"w") as file:
       file.write(objstr+"\n")
       file.write(consstr+"\n")
       file.write(boundstr+"\n")
       file.write("End\n")
    cplexoutpath = "{0}/{1}.lpout".format(runfolder,algoparstr)
    cplexscriptpath = "{0}/{1}.script".format(runfolder,algoparstr)
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
    print "QSAP History Inference done in {0} seconds".format(t2-t1)
    retvalues,objval = InputOutput.readCplexOutWithObj(cplexoutpath,specific = ["x"])
    assert HistoryUtil.checkValidCplexRun(cplexoutpath)
    for cplexfile in myutil.listfiles(runfolder):
        os.system("rm -rf {0}/{1}".format(runfolder,cplexfile))
    retkeys = retvalues.keys()    
    for key in retkeys:
        if retvalues[key] in [0.0,-0.0]:
           del retvalues[key]          
    return retvalues,objval


def runFlow(cutG,s,t,runfolder,sortedlabels,node2labels):
    """runs kolmogorov flow given cut graph
    Args:
       cutG:
       s:
       t:   
       runfolder:
       sortedlabels:
       node2labels:
    Returns:
       segment:  
       flow: 
    """
    assert type(cutG) == nx.DiGraph
    node2index,index2node = {},{}
    cutnodes = cutG.nodes()
    for index in xrange(len(cutnodes)):
        node = cutnodes[index]
        node2index[node] = index 
        index2node[index] = node
    mappedG = nx.DiGraph()
    for node in cutG.nodes():
        mappedG.add_node(node2index[node])
    for node1,node2 in cutG.edges():
        mappedG.add_edge(node2index[node1],node2index[node2],weight=cutG[node1][node2]["weight"])
    CODEPATH = "{0}/{1}".format(globals()["BOYKOVCODEDIR"],globals()["BOYKOVCODENAME"])
    mymap ={"s":s,"t":t}
    inputfile = "{0}/data{1}.input".format(runfolder,random.random())
    outfile = "{0}/data{1}.out".format(runfolder,random.random())
    cutUtil.storeBoykovGraph(mappedG,node2index[s],node2index[t],inputfile)
    code = "{0} {1} {2}".format(CODEPATH,inputfile,outfile)
    os.system(code)
    node2terminal,count,flow = {},0, None
    with open(outfile,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if count == 0:
               flow = float(line) 
               count +=1
               continue
            node,terminal = line.split("\t")
            realnode = index2node[int(node)]
            node2terminal[realnode] = mymap[terminal]
    node2assign = {}
    for node in node2labels.keys():
        foundlabel = None
        for index in xrange(len(sortedlabels)):
            label = sortedlabels[index]
            newnode = (node,label) 
            if node2terminal[newnode] == t:
               assert index != 0
               foundlabel = sortedlabels[index-1] 
               break
        if foundlabel == None:
           foundlabel = sortedlabels[-1]  
        node2assign[node] = foundlabel
        
    if globals()["isTest"]:
       assert flow <= 100000000000.0 
       QsapTest.checkFlowResult(mappedG,node2terminal,node2labels,sortedlabels,node2assign,s,t)
    for filename in [inputfile,outfile]:
        os.system("rm -rf {0}".format(filename))    
    return node2assign,flow
    

def addDataEdges(cutG,sortedlabels,node2labels,linear,s,t,INFWEIGHT):
    """adds data and constraint edges to cut graph
    Args:
      cutG: 
      sortedlabels:
      node2labels:
      linear:
      s:
      t:
      INFWEIGHT:
    Returns:
    """   
    mylinear = deepcopy(linear)
    minlinear = min([0.0] + linear.values())
    if minlinear < 0.0:
       for key in linear.keys():
           mylinear[key] += -1.0 * minlinear
    for node in node2labels.keys():
        curlabels = [s] + sortedlabels + [t]
        for index in xrange(1,len(curlabels)):
            prelabelstr,curlabelstr = curlabels[index-1:index+1]
            if prelabelstr == s:
               cutG.add_edge(s,(node,curlabelstr),weight=INFWEIGHT) 
            elif prelabelstr not in node2labels[node]:
               if curlabelstr == t:
                  cutG.add_edge((node,prelabelstr),t,weight=INFWEIGHT) 
               else:    
                  cutG.add_edge((node,prelabelstr),(node,curlabelstr),weight=INFWEIGHT)
            else: 
               varname = "x{0}?{1}".format(node,prelabelstr) 
               curval = 0.0
               if mylinear.has_key(varname):
                  curval = mylinear[varname]
               if curlabelstr == t:
                  cutG.add_edge((node,prelabelstr),t,weight=curval)
               else:    
                  cutG.add_edge((node,prelabelstr),(node,curlabelstr),weight=curval)
            if prelabelstr != s and curlabelstr != t:
               cutG.add_edge((node,curlabelstr),(node,prelabelstr),weight=INFWEIGHT)

               
def addInteractionEdges(G,cutG,node2labels,sortedlabels,quad):
    """adds interaction edges
    Args:
       G:
       cutG:
       node2labels:
       sortedlabels:
       quad:
    Returns:
    """ 
    undiredges,diredges = set(),set()
    for node1,node2 in G.edges():
        if G.has_edge(node2,node1) and (node2,node1) not in undiredges:
           undiredges.add((node1,node2))   
        elif not G.has_edge(node2,node1):
           diredges.add((node1,node2))
    assert len(diredges) + 2*len(undiredges) == G.number_of_edges()
    useedges = [("u",edge) for edge in undiredges] + [("d",edge) for edge in diredges]
    for sym,(node1,node2) in useedges:
        if not node2labels.has_key(node1) or not node2labels.has_key(node2):
           continue
        if len(node2labels[node1]) == 1 and len(node2labels[node2]) == 1:
           continue 
        
        mullen = len(sortedlabels)**2
        node1indices = [(index1,index2) for index1 in xrange(len(sortedlabels)) for index2 in xrange(1,len(sortedlabels))]
        node2indices = [(index1,index2) for index1 in xrange(len(sortedlabels)) for index2 in xrange(1,len(sortedlabels))]
        
        A = np.zeros((mullen,len(node1indices)+len(node2indices)),dtype=np.float64) 
        b = np.zeros((mullen,1),dtype=np.float64)
        delindices,keepindices = [],[]
        for index in xrange(mullen):
            in1 = index / len(sortedlabels)
            in2 = index % len(sortedlabels)
            label1 = sortedlabels[in1]
            label2 = sortedlabels[in2]
            if label1 not in node2labels[node1] or label2 not in node2labels[node2]:
               delindices.append(index)
               continue 
            keepindices.append(index)
            var1 = "x{0}?{1}".format(node1,label1)
            var2 = "x{0}?{1}".format(node2,label2)
            if sym == "d":
               if quad.has_key((var2,var1)):
                  b[index,0] = quad[(var2,var1)]
               else:   
                  b[index,0] = 0.0 
            elif sym == "u":
               if quad.has_key((var2,var1)):
                  b[index,0] = quad[(var2,var1)]
               elif quad.has_key((var1,var2)):
                  b[index,0] = quad[(var1,var2)]   
               else:   
                  b[index,0] = 0.0 
                     
            for curindex in xrange(len(node1indices)):
                index1,index2 = node1indices[curindex]
                if index1 <= in1 and index2 >= in2+1:
                   A[index,curindex] = 1.0 
            for curindex in xrange(len(node2indices)):
                index1,index2 = node2indices[curindex]  
                if index1 <= in2 and index2 >= in1+1:
                   A[index,len(node1indices)+curindex] = 1.0  
        b = [b[index,0] for index in xrange(np.shape(b)[0])]           
        assert len(delindices) == np.shape(A)[0] - (len(node2labels[node1]) *  len(node2labels[node2]))
        for colindex in xrange(len(node1indices)+len(node2indices)):
            assert A[mullen-1,colindex] == 0  
        if len(delindices) != 0:
           A = np.delete(A,delindices, 0)
           b = np.delete(b,delindices, 0)
           
        if True:
           print "seir info"
           print len(node2labels[node1])
           print len(node2labels[node2])
           print node2labels[node1]
           print node2labels[node2]
           #arank = np.linalg.matrix_rank(A)
           #print A
           #print b
           #print arank

        res = scipy.optimize.nnls(A,b)[0]
        if globals()["isTest"]:
           for index in xrange(np.shape(A)[0]):
               mysum = sum([A[index,index2]*res[index2] for index2 in xrange(np.shape(A)[1])])
               assert abs(mysum - b[index]) <= 0.0000001 
           neglist = [index for index in xrange(len(res)) if res[index] < -0.0001]
           assert len(neglist) == 0  
        
        for index in xrange(len(res)): 
            if abs(res[index]) <= math.exp(-11):
               continue 
            if index < len(node1indices):
               in1,in2 = node1indices[index] 
               sendnode,recnode = node1,node2
               sendlabel,reclabel = sortedlabels[in1],sortedlabels[in2]
            elif index >= len(node1indices):
               newindex = index - len(node1indices)
               in1,in2 = node2indices[newindex] 
               sendnode,recnode = node2,node1
               sendlabel,reclabel = sortedlabels[in1],sortedlabels[in2]
            assert not cutG.has_edge((sendnode,sendlabel),(recnode,reclabel))  
            cutG.add_edge((sendnode,sendlabel),(recnode,reclabel),weight=res[index]) 

               
def transform2Flow(G,node2vars,linear,quad,s,t,smodel,APPROX):
    """transforms graph into new graph for flow problem
    Args:
       G:
       node2vars:
       linear:
       quad: 
       s:
       t:
       smodel:
       APPROX:
    Returns:
       cutG:
    """
    INFWEIGHT = (sum(linear.values()) + sum([abs(val) for val in quad.values()])) * 1000 
    labelG,node2labels = getLabelDependency(node2vars,smodel,"weak")
    sortedlabels,sortednode2labels = sortLabels(labelG,node2labels)

    cutG = nx.DiGraph()
    cutG.add_node(s)
    cutG.add_node(t)
    for node in node2labels.keys(): 
        for labelstr in sortedlabels:
            cutG.add_node((node,labelstr))
    addDataEdges(cutG,sortedlabels,node2labels,linear,s,t,INFWEIGHT)
    if APPROX in ["reliability","prob"]:
       addInteractionEdges(G,cutG,node2labels,sortedlabels,quad)
       
    if globals()["isTest"]:  
       QsapTest.cutGraphCheck(cutG,sortedlabels,node2labels,G,quad,s,t,INFWEIGHT) 
    return cutG,sortedlabels,node2labels


def compareLabels(label1,label2,smodel,status):
    """Compares labels 
    Args:
       label1:
       label2:
       smodel: 
       status: weak or normal
    Returns:
       value: 1, 0, -1    
    """
    assert status in ["weak","normal"]
    if smodel in ["si","sir"]:
       sendstate,recstate = Trace.INFECTED, Trace.INFECTED
    elif smodel == "seir":
       sendstate,recstate = Trace.INFECTED, Trace.EXPOSED
    splitted1 = label1.split("?")
    times1 = {splitted1[2*index+0]: int(splitted1[2*index+1]) for index in xrange(len(splitted1)/2)}
    splitted2 = label2.split("?")
    times2 = {splitted2[2*index+0]: int(splitted2[2*index+1]) for index in xrange(len(splitted2)/2)}
    if HistoryUtil.checkAffect(times1,times2,smodel,sendstate,recstate,status):
       return 1 
    elif HistoryUtil.checkAffect(times2,times1,smodel,sendstate,recstate,status):
       return -1
    return 0  
             
def sortLabels(labelG,node2labels):
    """sort labels from given label dependency graph
    Args:
       labelG:
       node2labels:
    Returns:
       sortedlabels:
       sortednode2labels:
    """
    templabelG = nx.DiGraph(labelG)
    sortedlabels,sortednode2labels = [],{}        
    while templabelG.number_of_nodes() > 0:
        zeronode = None
        for node in templabelG.nodes():
            if len(templabelG.predecessors(node)) == 0:
               zeronode = node
               break  
        assert zeronode != None        
        sortedlabels.append(zeronode)
        templabelG.remove_node(zeronode)
    assert len(set(sortedlabels).symmetric_difference(set(labelG.nodes()))) == 0    
    labeldict = {sortedlabels[index]:index for index in xrange(len(sortedlabels))}
    for node in node2labels.keys():
        curlabels = list(node2labels[node])
        def label_cmp(label1, label2):
            return labeldict[label1] - labeldict[label2]
        sortednode2labels[node] = list(sorted(curlabels, cmp=label_cmp)) 
    if globals()["isTest"]:
       for index1 in xrange(len(sortedlabels)):
           node1 = sortedlabels[index1]
           for index2 in xrange(index1+1,len(sortedlabels)):
               node2 = sortedlabels[index2]
               assert not labelG.has_edge(node2,node1)      
    return sortedlabels,sortednode2labels


def getLabelDependency(node2vars,smodel,status):
    """returns label dependencies graph
    Args:
       node2vars:
       smodel:
       status: weak or normal
    Returns:
       labelG:
       node2labels:
    """ 
    assert status in ["normal","weak"]
    labelG = nx.DiGraph()
    node2labels = {}
    labels = set()
    for node in node2vars.keys():
        node2labels[node] = set()
        for var in node2vars[node]:
            labelstr = "?".join(var.split("?")[1:])
            labels.add(labelstr)
            node2labels[node].add(labelstr)
    for label in labels:
        labelG.add_node(label)
    print "info about labels"    
    print labels     
    for label1 in labels:
        for label2 in labels:
            if label1 == label2:
               continue
            comp = compareLabels(label1,label2,smodel,status)
            if comp == 1:
               labelG.add_edge(label1,label2)
            elif comp == -1:
               labelG.add_edge(label2,label1) 
    if globals()["isTest"]:
       assert len(nx.simple_cycles(labelG)) == 0 and nx.is_weakly_connected(labelG) 
       QsapTest.transDependencyCheck(labelG,status)
    return labelG,node2labels
    

def flow2Sol(node2assign,maxobstime):
    """converts flow results to sol
    Args:
       node2assign:
       maxobstime:
    Returns:
       sol:   
    """
    states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    sol = {time: {state:set() for state in states} for time in xrange(maxobstime+1)}
    for node in node2assign.keys():
        label = node2assign[node]
        splitted = label.split("?")
        times = {splitted[2*index+0]: int(splitted[2*index+1]) for index in xrange(len(splitted)/2)}
        for state in times.keys():
            sol[times[state]][state].add(node)
    return sol


def assignVar2label(node2var):
    """converts var to label assignment
    Args: 
       node2var:
    Returns:
       node2label:
    """
    node2label = {}
    for node in node2var.keys():
        label = node2var[node].replace("x{0}?".format(node),"")
        node2label[node] = label
    return node2label
    

def estimateRatios(retvalues,node2vars,linear,quad,objval):
    """estimates rounding ratio
    Args:
       retvalues:
       node2vars:
       linear:
       quad:
       objval:
    Returns:
       minratio,maxratio,avgratio:
    """
    singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
    doublevars = [var for var in retvalues.keys() if var.find("_") != -1]
    vals = set(retvalues.values())
    print "seenvals",vals
    valmap = {val:0 for val in vals}
    for singlevar in singlevars:
        valmap[retvalues[singlevar]] += 1
    print valmap    
    roundvals = []
    loopcount = 100
    for index in xrange(loopcount):
        node2assign = {}
        for node in node2vars.keys():
            curvars = list(node2vars[node])
            random.shuffle(curvars)
            randval = random.random()
            cursum = 0.0
            foundvar = None
            for var in curvars:
                if retvalues.has_key(var):
                   cursum += retvalues[var]
                if cursum >= randval:
                   foundvar = var
                   break
            if foundvar == None:
               foundvar = curvars[-1]        
            label = foundvar.replace("x{0}?".format(node),"")  
            node2assign[node] = label
        localval = estimateObjVal(node2assign,linear,quad,G)
        roundvals.append(localval)
    minratio,maxratio,avgratio = objval/min(roundvals), objval/max(roundvals), objval/(sum(roundvals)/loopcount)
    return minratio,maxratio,avgratio

def randomRound2Sol(retvalues,node2vars,maxobstime):
    """does randomized rounding
    Args:
       retvalues:
       node2vars:
       maxobstime:
    Returns:
       sol:
       node2var:
    """
    states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    sol = {time: {state:set() for state in states} for time in xrange(maxobstime+1)}            
    node2var = {}
    for node in node2vars.keys():
        curvars = list(node2vars[node])
        random.shuffle(curvars)
        randval = random.random()
        cursum = 0.0
        foundvar = None
        for var in curvars:
            if retvalues.has_key(var):
               cursum += retvalues[var]
            if cursum >= randval:
               foundvar = var
               break
        if foundvar == None:
           foundvar = curvars[-1]
        node2var[node] = foundvar   
        nodetimes = HistoryUtil.var2Times(foundvar)
        for state in nodetimes.keys():
            sol[nodetimes[state]][state].add(node)
    return sol,node2var
 
FLEXIBLE = True
def runLp(modcurstates,smodel,G,runfolder,node2slots,node2vars,linear,quad,ALGOTYPE,algoinfo):
    """runs lp (randomized rounds for frac values of ALGOTYPE=special type)
    Args:
       modcurstates:
       smodel:
       G:
       runfolder:
       node2slots:
       node2vars:
       linear:
       quad:
       ALGOTYPE: general or special
       algoinfo:
    Returns:
       sol:
       retobjval: return objval
    """
    APPROX = algoinfo["approx"]
    if algoinfo["approx"] == "ratio":
       assert algoinfo.has_key("roundmethod") 
    maxtime = max(modcurstates.keys())
    objstr = genMLObj(linear,quad,algoinfo)
    consstr = genMLMainConst(node2vars,linear,smodel,G,quad)
    boundstr = genMLBoundConst(linear,quad,node2vars,G)
    retvalues,objval,runtime = HistoryUtil.runHistoryCode(consstr,objstr,boundstr,runfolder,["x"],"qsap-{0}".format(APPROX),"removeZeros")
     
    if APPROX in ["reliability","prob"] and ALGOTYPE == "special":
       isfrac = False 
       singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
       for singlevar in singlevars:
           if retvalues[singlevar] != 1.0:
              isfrac = True
              break
       if not globals()["FLEXIBLE"]:
          assert not isfrac         
    else:
       isfrac = True  
    if isfrac:
       if APPROX == "ratio" and algoinfo["roundmethod"] == "random":
          sol,node2var = randomRound2Sol(retvalues,node2vars,maxtime)
       else:
          sol,node2var = randomRound2Sol(retvalues,node2vars,maxtime)    
       node2label = assignVar2label(node2var)  
       retobjval = estimateObjVal(node2label,linear,quad,G)
       #print "ratios"
       #print objval,retobjval
       #assert objval != 0.0 and retobjval != 0.0
       #ratio = objval /retobjval
       #print ratio
       #assert ratio <= 1.5
       #exit(1)
    else:
       sol = convert2Sol(retvalues,maxtime)
       retobjval = objval
       node2label = {}
       for singlevar in singlevars:
           node = int(singlevar.replace("x","").split("?")[0])
           label = singlevar.replace("x{0}?".format(node),"")
           node2label[node] = label
                    
    if globals()["isTest"]:
       paraminfo = [node2slots,node2vars,modcurstates,quad,linear,G,smodel]
       solinfo = [retvalues,sol,objval,retobjval,isfrac,node2label] 
       QsapTest.lpSolutionCheck(paraminfo,solinfo)
    return sol,node2label,retobjval,objval 
     

RUNMODE = "lp"  #"flow","lp","test"
BOYKOVCODEDIR = "maxflow-v3.02"
BOYKOVCODENAME = "boykovcut"
def runQsap(modcurstates,smodel,G,runfolder,algoinfo,ALGOTYPE):
    """Qsap relaxation for known interval
       Args:
         modcurstates: current state knowledge
         smodel:
         G: Graph
         runfolder:
         algoinfo:
         ALGOTYPE:
       Returns:
         sol: solution
    """
    APPROX = algoinfo["approx"]
    maxtime = max(modcurstates.keys())
    node2slots = assignStateTimeSlots(modcurstates,smodel,G)
    node2vars = getNode2Vars(node2slots,smodel)
    print modcurstates.keys()
    for time in sorted(modcurstates.keys()):
        print time
        print modcurstates[time][Trace.SUSCEPTIBLE]
        print modcurstates[time][Trace.EXPOSED]
        print modcurstates[time][Trace.INFECTED]
        print modcurstates[time][Trace.RECOVERED]
    linear,quad = estimateCoefs(node2vars,smodel,G,maxtime,APPROX)  
    #if APPROX == "ratoo"
    if ALGOTYPE == "special":
       if globals()["RUNMODE"] in ["flow","test"]:
          s,t = -1, max(G.nodes()) + 1
          cutG,sortedlabels,node2labels = transform2Flow(G,node2vars,linear,quad,s,t,smodel,APPROX)
          node2assign,flow = runFlow(cutG,s,t,runfolder,sortedlabels,node2labels)
          print "flow is;",flow
          assert len(node2assign.keys()) == len(node2vars.keys())
          flowsol = flow2Sol(node2assign,maxtime)
       if globals()["RUNMODE"] in ["lp","test"]:
          lpsol,lpnode2assign,retobjval,objval  = runLp(modcurstates,smodel,G,runfolder,node2slots,node2vars,linear,quad,ALGOTYPE,algoinfo)
       if globals()["RUNMODE"] in ["test"]:
          print retobjval,objval
          objlp = estimateObjVal(lpnode2assign,linear,quad,G)
          objflow = estimateObjVal(node2assign,linear,quad,G)
          print objlp,objflow,flow
          assert len(node2assign.keys()) == len(lpnode2assign.keys())
          assert objflow <= objlp and abs(objval - objflow ) <= 0.01 and ( abs(objflow - flow) <= 0.01 or flow <= objflow) #flow ignores the interaction when there is only single variable!!
          for index in xrange(10):
              node2assign = {}
              for node in node2vars.keys():
                  curvars = list(node2vars[node])
                  random.shuffle(curvars)
                  varname = curvars[0]
                  label = varname.replace("x{0}?".format(node),"")
                  node2assign[node] = label
              curobjval = estimateObjVal(node2assign,linear,quad,G)
              assert curobjval >= objlp and curobjval >= objflow
    elif ALGOTYPE == "general": 
       if APPROX == "ratio":
          pass             
       elif APPROX in ["reliability","prob"]:
          lpsol,lpnode2assign,retobjval,objval  = runLp(modcurstates,smodel,G,runfolder,node2slots,node2vars,linear,quad,ALGOTYPE,algoinfo)
           
    if globals()["isTest"]:
       [QsapTest.checkVariableValidity(modcurstates,varname,smodel) for node in node2vars.keys() for varname in node2vars[node]]
       QsapTest.coefCheck(quad,linear,smodel,G) 
    if ALGOTYPE == "general" or globals()["RUNMODE"] == "lp":
       return lpsol
    elif globals()["RUNMODE"] in ["flow","test"]:
       return flowsol  

    

def sortVars(varnames):
    """sorts vars according to time
    Args:
       varnames:
    Returns:
       sortednames:
    """
    esorteddict,isorteddict = {},{}
    for varname in varnames:
        ttimes = HistoryUtil.var2Times(varname)
        if ttimes.has_key(Trace.INFECTED):
           isorteddict.setdefault(ttimes[Trace.INFECTED],set()).add(varname)
        elif ttimes.has_key(Trace.EXPOSED):
           esorteddict.setdefault(ttimes[Trace.EXPOSED],set()).add(varname)   
    sortedvars = []
    for time in sorted(isorteddict.keys()):
        sortedvars.extend(isorteddict[time]) 
    for time in sorted(esorteddict.keys()):
        sortedvars.extend(esorteddict[time])     
    return sortedvars


def convert2Sol(retvalues,maxobstime):
    """converts retvalues to sol format
    Args:
       retvalues:
       maxobstime:
    Returns:
       sol:
    """
    print "passed from here!!"
    singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
    states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    sol = {time: {state:set() for state in states} for time in xrange(maxobstime+1)}
    for varname in singlevars:
        node = int(varname.replace("x","").split("?")[0])
        nodetimes = HistoryUtil.var2Times(varname)
        for state in nodetimes.keys():
            sol[nodetimes[state]][state].add(node)
    print "solinfo:"        
    for time in sol.keys():
        print time
        print sol[time][Trace.SUSCEPTIBLE]
        print sol[time][Trace.EXPOSED]
        print sol[time][Trace.INFECTED]
        print sol[time][Trace.RECOVERED]      
    return sol
    
def getAlgoType(G,smodel,APPROX):
    """gets algo type
    Args:
       G:
       smodel:
       APPROX:
    Returns:
       ALGOTYPE: 
    """    
    ALGOTYPE = "special"
    if smodel in ["si","sir"]:
       transstr = Trace.S2I
    elif smodel == "seir":
       transstr = Trace.S2E 
    if APPROX == "ratio":
       ALGOTYPE = "general"
    elif APPROX in ["reliability","prob"]:     
       for node1,node2 in G.edges():
           if G[node1][node2][transstr][0] not in ["expo","rayleigh","normal","geometric"]:
              ALGOTYPE = "general"
              break
    return ALGOTYPE  


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
    globals()["isTest"] = TESTMODE
    assert checkUtil.checkPreTrace(curstates,obstimes)
    seenstates = deepcopy(curstates)
    print "algoinfo"
    print algoinfo
    exit(1)
    assert algoinfo["approx"] in ["prob","reliability","ratio"]
    ensemble = algoinfo["ensemble"]
    APPROX = algoinfo["approx"]
    enscount = 1
    if ensemble:
       enscount = algoinfo["enscount"]
    allsol = []
    ALGOTYPE = getAlgoType(G,smodel,APPROX)
    for index in xrange(enscount):
        modcurstates = {obstimes[timeindex]:deepcopy(seenstates[timeindex]) for timeindex in xrange(len(obstimes))}
        history = runQsap(modcurstates,smodel,G,runfolder,algoinfo,ALGOTYPE)
        history = HistoryUtil.lastFixSolution(history,G,smodel,"bound",infermode) 
        allsol.append(deepcopy(history)) 
    if ensemble:
       history = DisEns.getEnsembleSolution(allsol,infermode,seenstates,G,smodel)
    if infermode == "Spreader":
       foundtime = None
       for time in sorted(allsol[0].keys()):
           if len(allsol[0][time][Trace.INFECTED]) != 0:
              foundtime = time
              break 
       history = allsol[0][foundtime], 0
       print "Inferred Spreader Solution at time {0}".format(history[1])
       print history[0]
    elif infermode == "History":     
       history = allsol[0]
       print "modcurstates:"
       for time in sorted(modcurstates.keys()):  
           print time,len(modcurstates[time][Trace.SUSCEPTIBLE]),len(modcurstates[time][Trace.EXPOSED]),len(modcurstates[time][Trace.INFECTED]),len(modcurstates[time][Trace.RECOVERED])        
       print "inferred history solution"
       for time in sorted(history.keys()):
           print time,len(history[time][Trace.SUSCEPTIBLE]),len(history[time][Trace.EXPOSED]),len(history[time][Trace.INFECTED]),len(history[time][Trace.RECOVERED])
       
    if globals()["isTest"]:
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
