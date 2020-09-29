"""Bounded History Reconstruction Utilities
"""
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


def sortLabels(labelG,node2labels):
    """sort labels from given label dependency graph
    Args:
       labelG:
       node2labels:
    Returns:
       sortedlabels:
       sortednode2labels:
    """
    TESTMODE = True
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
    if TESTMODE:
       for index1 in xrange(len(sortedlabels)):
           node1 = sortedlabels[index1]
           for index2 in xrange(index1+1,len(sortedlabels)):
               node2 = sortedlabels[index2]
               assert not labelG.has_edge(node2,node1)      
    return sortedlabels,sortednode2labels

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
    TESTMODE = True
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
    for label1 in labels:
        for label2 in labels:
            if label1 == label2:
               continue
            comp = compareLabels(label1,label2,smodel,status)
            if comp == 1:
               labelG.add_edge(label1,label2)
            elif comp == -1:
               labelG.add_edge(label2,label1) 
    if TESTMODE:
       assert len(nx.simple_cycles(labelG)) == 0 and nx.is_weakly_connected(labelG) 
       transDependencyCheck(labelG,status)
    return labelG,node2labels


def checkValidityDirectly(sol,G,labelgraph,mode):
    """checks valid diffusion validity directly not over labelnode graph
    Args:
       sol:
       G:
       labelgraph: dependencies between labels(not nodes)
       mode: var/label
    Returns:
       bool: true/false
    """
    assert mode in ["label","var"]
    labelsol = deepcopy(sol)
    if mode == "var":
       for node in labelsol.keys(): 
           label = labelsol[node].replace("x{0}?".format(node),"")
           labelsol[node] = label
    for node in labelsol.keys():
        label = labelsol[node]
        if len(labelgraph.predecessors(label)) == 0:
           continue    
        flag = False 
        for givernode in set(labelsol.keys()).intersection(G.predecessors(node)):
            if labelgraph.has_edge(labelsol[givernode],label):
               flag = True
               break
        if not flag:
           return False      
    return True
    

def assignLabel2Var(node2label):
    """converts var to label assignment
    Args: 
       node2label:
    Returns:
       node2var:
    """
    node2var = {}
    for node in node2label.keys():
        var = "x{0}?{1}".format(node,node2label[node])
        node2var[node] = var
    return node2var

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


def makeLabelNodeGraph(G,node2vars,smodel,comparemode,nodemode="label"):
    """makes big label-node or var-node graph to use in dependency checking
    Args:
       G:
       node2vars:
       smodel:
       comparemode:
       nodemode:
    Returns:
       labelnodeG:
    """
    assert nodemode in ["label","var"]
    assert comparemode in ["normal","weak"]
    labelnodeG = nx.DiGraph()
    for node in node2vars.keys():
        for var in node2vars[node]:
            sym = var
            if nodemode == "label":
               sym = var.replace("x{0}?".format(node),"")
            labelnodeG.add_node((node,sym))
    for node1,sym1 in labelnodeG.nodes():
        label1 = sym1
        if nodemode == "var":
           label1 = sym1.replace("x{0}?".format(node1),"")
        for node2,sym2 in labelnodeG.nodes():
            if node1 == node2:
               continue
            label2 = sym2 
            if nodemode == "var":
               label2 = sym2.replace("x{0}?".format(node2),"")
            if G.has_edge(node1,node2) and compareLabels(label1,label2,smodel,comparemode) == 1:
               labelnodeG.add_edge((node1,sym1),(node2,sym2)) 
    return labelnodeG



def estimateTransCoefs(node,times1,varname1,smodel,interinfo,G,node2vars,ZERO,APPROX,affectmode):
    """estimate transition coefs(s2i or s2e)
    Args:
       node:
       times1:
       varname1:
       smodel:
       interinfo:
       G:
       node2vars:
       ZERO:
       APPROX:
       affectmode:
    Returns:
       quad:
       linear:
    """
    sstate,rstate,transdist = interinfo 
    assert affectmode in ["weak","normal"]
    quad,linear = {}, {} 
    for neighnode in set(G.predecessors(node)).intersection(set(node2vars.keys())):
        for varname2 in node2vars[neighnode]:
            times2 = HistoryUtil.var2Times(varname2)
            if HistoryUtil.checkAffect(times2,times1,smodel,sstate,rstate,affectmode):
               f1 = Dist.genPartDist(G[neighnode][node][transdist][0],G[neighnode][node][transdist][1],"reverseCdf",G[neighnode][node][Trace.SPROB])
               fdist = Dist.genPdf(f1)
               fdist[0] = 1.0
               p1 = Dist.genPartDist(G[neighnode][node][transdist][0],G[neighnode][node][transdist][1],"normal",G[neighnode][node][Trace.SPROB])
               pdist = Dist.genPdf(p1)
               diftime = times1[rstate] - times2[sstate]
               assert diftime > 0
               if APPROX == "taylor":
                  part1,part2 = ZERO,ZERO
                  if fdist.has_key(diftime) and fdist[diftime] > ZERO:
                     part1 = fdist[diftime]
                  if fdist.has_key(diftime-1) and fdist[diftime-1] > ZERO:
                     part2 = fdist[diftime-1]
                  val = part2**2 /part1
               elif APPROX == "reliability": 
                  val = ZERO
                  if fdist.has_key(diftime-1) and fdist[diftime-1] > ZERO:
                     val = fdist[diftime-1] 
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
               if val != 1.0:  #if abs(val-1.0) >= 0.000001:
                  quad[(varname1,varname2)] = -1.0*math.log(val)       
               #if APPROX == "taylor" and linearval != 1.0:
               #   assert linearval <= 1.0
               #   linear.setdefault(varname1,0.0)
               #   linear[varname1] += math.log(linearval)              
    return quad,linear


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


def estimateCoefs(node2vars,smodel,G,maxtime,APPROX,affectmode):
    """estimate coefficients
    Args:
      node2vars:
      smodel:
      G:
      maxtime:
      APPROX:
      affectmode:
    Returns:
      linear:
      quad: 
    """
    assert APPROX in ["taylor","prob","ratio","reliability"]  
    ZERO = 0.000000000000000001 
    if smodel in ["si","sir"]:
       interinfo = Trace.INFECTED,Trace.INFECTED,Trace.S2I
    elif smodel == "seir":
       interinfo = Trace.INFECTED,Trace.EXPOSED,Trace.S2E
    linear = {varname:0.0 for node in node2vars.keys() for varname in node2vars[node]}
    quad = {} 
    snodes = set(G.nodes()).difference(set(node2vars.keys()))
    for node in node2vars.keys():
        for varname1 in node2vars[node]:
            times1 = HistoryUtil.var2Times(varname1)
            localquad,locallinear = estimateTransCoefs(node,times1,varname1,smodel,interinfo,G,node2vars,ZERO,APPROX,affectmode)
            for curvarname1,varname2 in localquad.keys(): 
                assert not quad.has_key((varname1,varname2)) 
                quad[(varname1,varname2)] = localquad[(varname1,varname2)]
            for curvarname in locallinear.keys():
                linear[curvarname] += locallinear[curvarname]     
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
               for sidevar in set(node2vars[node]).difference(set([mainvar])):
                   applinear.setdefault(sidevar,0.0)
                   applinear[sidevar] += linear[mainvar]
       linear = deepcopy(applinear)  
    return linear,quad

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
    TESTMODE = True
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
    
    if TESTMODE:     
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


def getNode2Vars(node2slots,smodel):
    """gets node2vars
    Args:
       node2slots:
       smodel:
    Returns:
       node2vars:
    """
    TESTMODE = True
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
                  
    if TESTMODE:
       checkNodeVars(node2vars,smodel)
    return node2vars

