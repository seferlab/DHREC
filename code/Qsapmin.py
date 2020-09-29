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
import boundInferUtil

class QsapMinTest():
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
     def coefCheck(quad,linear,smodel,G,affectmode):
         """coef check 
         Args:
           quad:
           linear: 
           smodel:
           G:
           affectmode:
         Returns:
           bool: true/false
         """
         assert affectmode in ["weak","normal"]
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
             assert HistoryUtil.checkAffect(times2,times1,smodel,sstate,rstate,affectmode) and G.has_edge(node2,node1) 
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
         INFTY = 100000000
         if smodel in ["si","sir"]:
            if not nodetimes.has_key(Trace.RECOVERED):
               nodetimes[Trace.RECOVERED] = INFTY
            for time in sortedtimes:
                if time < nodetimes[Trace.INFECTED]:
                   assert node in modcurstates[time][Trace.SUSCEPTIBLE]
                elif time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                   assert node in modcurstates[time][Trace.INFECTED]
                elif time >= nodetimes[Trace.RECOVERED]:
                   assert node in modcurstates[time][Trace.RECOVERED]
         elif smodel == "seir":
            if not nodetimes.has_key(Trace.EXPOSED):
               mode = "imode"
               assert nodetimes[Trace.INFECTED] == 0 
               if not nodetimes.has_key(Trace.RECOVERED):
                  nodetimes[Trace.RECOVERED] = INFTY
            else:
               mode = "emode"
               assert nodetimes[Trace.EXPOSED] > 0 
               if not nodetimes.has_key(Trace.INFECTED):
                  nodetimes[Trace.INFECTED] = INFTY
               if not nodetimes.has_key(Trace.RECOVERED):
                  nodetimes[Trace.RECOVERED] = INFTY
            if mode == "emode":   
               for time in sortedtimes:
                   if time < nodetimes[Trace.EXPOSED]:
                      assert node in modcurstates[time][Trace.SUSCEPTIBLE]
                   elif time >= nodetimes[Trace.EXPOSED] and time < nodetimes[Trace.INFECTED]:
                      assert node in modcurstates[time][Trace.EXPOSED]
                   elif time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.INFECTED]   
                   elif time >= nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.RECOVERED]
            elif mode == "imode":
               for time in sortedtimes:
                   if time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.INFECTED]
                   elif time >= nodetimes[Trace.RECOVERED]:
                      assert node in modcurstates[time][Trace.RECOVERED]       
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

     
def estimateObjVal(node2assign,linear,quad,G):
    """estimate objective value
    Args:
      node2assign:
      linear:
      quad:
      G:
    Returns:
      objval: 
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
    for node in node2vars.keys():
        for var1 in node2vars[node]:
            tnodes = set(G.predecessors(node)).union(set(G.successors(node)))
            affectnodes = tnodes.intersection(node2vars.keys())
            for neighnode in affectnodes:
                for var2 in node2vars[neighnode]:
                    boundstr += "0 <= {0}_{1} <= 1.0\n".format(var1,var2)
                    boundstr += "0 <= {0}_{1} <= 1.0\n".format(var2,var1)
    return boundstr

def genMLObj(linear,quad):
    """generates QSAP objective function 
    Args:
      linear:
      quad:
    Returns:
      objstr:
    """       
    isPlus = lambda x: "+" if x >= 0 else " "
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
            affectnodes = tnodes.intersection(node2vars.keys()) #affectnodes = set(G.predecessors(node)).intersection(node2vars.keys())
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
       QsapMinTest.checkFlowResult(mappedG,node2terminal,node2labels,sortedlabels,node2assign,s,t)
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

               
def addInteractionEdges2(G,cutG,node2labels,sortedlabels,quad):
    """adds interaction edges
    Args:
       G:
       cutG:
       node2labels:
       sortedlabels:
       quad:
    Returns:
    """
    SUMVAL = 0.0
    if len(quad.values()) > 0 and min(quad.values()) < 0.0:
       SUMVAL = -1.0*min(quad.values())
    undiredges,diredges = set(),set()
    for node1,node2 in G.edges():
        if G.has_edge(node2,node1) and (node2,node1) not in undiredges:
           undiredges.add((node1,node2))   
        elif not G.has_edge(node2,node1):
           diredges.add((node1,node2))
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
        b = [b[index,0]+SUMVAL for index in xrange(np.shape(b)[0])]
        #print b
        assert len(delindices) == np.shape(A)[0] - (len(node2labels[node1]) *  len(node2labels[node2]))
        for colindex in xrange(len(node1indices)+len(node2indices)):
            assert A[mullen-1,colindex] == 0  
        if len(delindices) != 0:
           A = np.delete(A,delindices, 0)
           b = np.delete(b,delindices, 0)
           
        if False:
           print " info"
           print len(node2labels[node1])
           print len(node2labels[node2])
           print node2labels[node1]
           print node2labels[node2]
           print np.linalg.matrix_rank(A)
           print np.shape(A)
           print np.shape(b)
        
        res = scipy.optimize.nnls(A,b)[0]
        #print res
        if globals()["isTest"]:
           for index in xrange(np.shape(A)[0]):
               mysum = sum([A[index,index2]*res[index2] for index2 in xrange(np.shape(A)[1])])
               #print mysum,b[index],abs(mysum-b[index])
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
    SUMVAL = 0.0
    if len(quad.values()) > 0 and min(quad.values()) < 0.0:
       SUMVAL = -1.0*min(quad.values())
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
        b = [b[index,0]+SUMVAL for index in xrange(np.shape(b)[0])]
        #print b
        assert len(delindices) == np.shape(A)[0] - (len(node2labels[node1]) *  len(node2labels[node2]))
        for colindex in xrange(len(node1indices)+len(node2indices)):
            assert A[mullen-1,colindex] == 0  
        if len(delindices) != 0:
           A = np.delete(A,delindices, 0)
           b = np.delete(b,delindices, 0)
           
        if False:
           print " info"
           print len(node2labels[node1])
           print len(node2labels[node2])
           print node2labels[node1]
           print node2labels[node2]
           print np.linalg.matrix_rank(A)
           print np.shape(A)
           print np.shape(b)
        
        res = scipy.optimize.nnls(A,b)[0]
        #print res
        if globals()["isTest"]:
           for index in xrange(np.shape(A)[0]):
               mysum = sum([A[index,index2]*res[index2] for index2 in xrange(np.shape(A)[1])])
               #print mysum,b[index],abs(mysum-b[index])
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
    labelG,node2labels = boundInferUtil.getLabelDependency(node2vars,smodel,"weak")
    sortedlabels,sortednode2labels = boundInferUtil.sortLabels(labelG,node2labels)
    
    cutG = nx.DiGraph()
    cutG.add_node(s)
    cutG.add_node(t)
    for node in node2labels.keys(): 
        for labelstr in sortedlabels:
            cutG.add_node((node,labelstr))
    addDataEdges(cutG,sortedlabels,node2labels,linear,s,t,INFWEIGHT)
    if APPROX in ["reliability","prob","taylor"]:
       addInteractionEdges(G,cutG,node2labels,sortedlabels,quad)
       
    if globals()["isTest"]:  
       QsapMinTest.cutGraphCheck(cutG,sortedlabels,node2labels,G,quad,s,t,INFWEIGHT) 
    return cutG,sortedlabels,node2labels

            
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
    

def convert2Sol(retvalues,maxobstime):
    """converts retvalues to sol format
    Args:
       retvalues:
       maxobstime:
    Returns:
       sol:
    """
    singlevars = [var for var in retvalues.keys() if var.find("_") == -1]
    states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
    sol = {time: {state:set() for state in states} for time in xrange(maxobstime+1)}
    for varname in singlevars:
        node = int(varname.replace("x","").split("?")[0])
        nodetimes = HistoryUtil.var2Times(varname)
        for state in nodetimes.keys():
            sol[nodetimes[state]][state].add(node)
    return sol

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
    maxtime = max(modcurstates.keys())
    objstr = genMLObj(linear,quad)
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
       if ALGOTYPE == "special": 
          sol,node2var = randomRound2Sol(retvalues,node2vars,maxtime)    
          node2label = boundInferUtil.assignVar2label(node2var)  
          retobjval = estimateObjVal(node2label,linear,quad,G)
       elif ALGOTYPE == "general":
          print "Rounding for general case not implemented!!"
          exit(1)    
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
       QsapMinTest.lpSolutionCheck(paraminfo,solinfo)
    return sol,node2label,retobjval,objval 


RUNMODE = "flow"  #"flow","lp","test"
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
    node2slots = boundInferUtil.assignStateTimeSlots(modcurstates,smodel,G)
    node2vars = boundInferUtil.getNode2Vars(node2slots,smodel)
    linear,quad = boundInferUtil.estimateCoefs(node2vars,smodel,G,maxtime,APPROX,"weak") 
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
       lpsol,lpnode2assign,retobjval,objval  = runLp(modcurstates,smodel,G,runfolder,node2slots,node2vars,linear,quad,ALGOTYPE,algoinfo)
           
    if globals()["isTest"]:
       [QsapMinTest.checkVariableValidity(modcurstates,varname,smodel) for node in node2vars.keys() for varname in node2vars[node]]
       QsapMinTest.coefCheck(quad,linear,smodel,G,"weak") 
    if ALGOTYPE == "general" or globals()["RUNMODE"] == "lp":
       return lpsol
    elif globals()["RUNMODE"] in ["flow","test"]:
       return flowsol  
    
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
    transmap = {"si":Trace.S2I,"sir":Trace.S2I,"seir":Trace.S2E}
    if APPROX in ["reliability","prob"]:     
       for node1,node2 in G.edges():
           if G[node1][node2][transmap[smodel]][0] not in ["expo","rayleigh","normal","geometric"]:
              ALGOTYPE = "general"
              break
    return ALGOTYPE  


def checkParamSuitable(algoinfo,G,smodel):
    """checks whether algo is suitable or not
    Args:
       algoinfo:
       G:
       smodel:
    Returns:
       bool: true/false
    """
    transmap = {"si":Trace.S2I,"sir":Trace.S2I,"seir":Trace.S2E}
    if algoinfo["approx"] == "taylor":
       for node1,node2 in G.edges():
           if G[node1][node2][transmap[smodel]][0] not in ["expo","rayleigh","normal","geometric"]:
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
    globals()["isTest"] = TESTMODE
    assert algoinfo["approx"] in ["taylor","prob","reliability"] and checkParamSuitable(algoinfo,G,smodel) and checkUtil.checkPreTrace(curstates,obstimes)
    seenstates = deepcopy(curstates)
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
       history = DisEns.getEnsembleSolution(allsol,infermode,seenstates,G,smodel) #?
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
    exit(1)   
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
