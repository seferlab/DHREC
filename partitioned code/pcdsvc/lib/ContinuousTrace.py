import networkx as nx
import numpy as np
import scipy as sp
import random
import math
import sys
import os
import myutilities as myutil
import gzip
import string
from ContinuousDistribution import ContinuousDistribution as Dist
from Trace import Trace

class ContinuousTrace(Trace):
    """ Continuous Trace Generation class -> inherits from Trace
    """
    def __init__():
       return

    @classmethod
    def genTraces(self,G,startnodecount,evol,tracecount,smodel,modelparams={},startnode=None):
       """generates continous traces
       Args:
          G: Graph with attributes
          startnodecount: start node count
          evol: static/dynamic
          tracecount: number of traces
          smodel: spreading model
          modelparams: model parameters
          startnode: startnode
       Returns:
          traces: generated traces 
       """
       traces = []
       assert smodel in ["si","sir","seir","sis"]
       assert evol == "static"
       for index in xrange(tracecount):
           funcname = "genCont"
           method = getattr(ContinuousTrace,funcname)
           (stimes,etimes,itimes,rtimes) = method(G,startnodecount,smodel,evol,modelparams,startnode)
           trace = self.convertTimes2Trace(smodel,stimes,etimes,itimes,rtimes)
           traces.append(trace)
       return traces
   
    @classmethod
    def genCont(self,G,startnodecount,smodel,evol,modelparams,startnode=None):
       """generates traces from continous distribution
       Args:
         G: Graph with attributes on it
         startnodecount: start node count
         smodel: spreading model
         evol: static/dynamic graph
         modelparams: spreading model parameters(such as maxtime for sis)
         startnode: start node of trace
       Returns:
         trace: generated trace
       """
       allnodes = G.nodes()
       if startnode == None: 
          random.shuffle(allnodes)
          startnodes = allnodes[startnodecount:startnodecount+1]
       else:
          startnodes = [int(startnode)]
       stimes,etimes,itimes,rtimes = {}, {}, {}, {}
       iactive,eactive = {}, {}
       if smodel in ["si","sis","sir"]:
          iactive = {node: 0.0 for node in startnodes} #active in terms of infected state
       elif smodel == "seir":
          eactive = {node: 0.0 for node in startnodes}  #active in terms of exposed (used in seir)
       if smodel == "sis":
          self.MAXSPREADTIME = modelparams["sismaxtime"]
          stimes = {node: set([0.0]) for node in G.nodes() if node not in startnodes}
       while len(iactive.keys())!=0 or len(eactive.keys())!=0:
          minnode, mintime = None, 100000000000000000.0
          if len(iactive.keys()) != 0:     
             tnode, ttime = min(iactive.items(), key=lambda x: x[1])
             if ttime < mintime:
                minnode, mintime = tnode, ttime
                ie = self.INFECTED
          if smodel == "seir":
             if len(eactive.keys())!=0:  
                tnode, ttime = min(eactive.items(), key=lambda x: x[1])
                if ttime < mintime:
                   minnode, mintime = tnode, ttime
                   ie = self.EXPOSED
          if mintime > self.MAXSPREADTIME:
             break
          if ie == self.EXPOSED:
             edur = Dist.genContRandNum(G.node[minnode][self.E2I][0],G.node[minnode][self.E2I][1])
             assert not iactive.has_key(minnode)
             etimes[minnode] = mintime
             iactive[minnode] = mintime + edur
             del eactive[minnode]
             continue
          if smodel in ["si","sir","seir"]:
             itimes[minnode] = mintime
          elif smodel == "sis":
             itimes.setdefault(minnode,set())
             itimes[minnode].add(mintime)
             idur = Dist.genContRandNum(G.node[minnode][self.I2S][0],G.node[minnode][self.I2S][1])
             stimes.setdefault(minnode,set())
             stimes[minnode].add(mintime+idur)
          curtime = mintime
          del iactive[minnode]
          if smodel in ["sir","seir"]:
             rtimes[minnode] = curtime + Dist.genContRandNum(G.node[minnode][self.I2R][0],G.node[minnode][self.I2R][1])
          for node in G[minnode].keys():
              if minnode == node:
                 continue 
              if smodel in ["si","sir"] and itimes.has_key(node):
                 continue
              elif smodel == "seir" and (etimes.has_key(node) or itimes.has_key(node)):
                 continue
              elif smodel == "sis": #determine whether node infected or susceptible for sis model
                 imax, smax = -1, -1
                 if itimes.has_key(node):
                    imax = max(itimes[node])
                 if stimes.has_key(node):
                    smax = max(stimes[node])
                 if imax != -1 and smax != -1:
                    if smax > curtime and imax < curtime: #currently infected
                       continue
              if smodel in ["si","sir","sis"]:
                 affectlen = Dist.genContRandNum(G[minnode][node][self.S2I][0],G[minnode][node][self.S2I][1],G[minnode][node][self.SPROB])
              elif smodel == "seir":
                 affectlen = Dist.genContRandNum(G[minnode][node][self.S2E][0],G[minnode][node][self.S2E][1],G[minnode][node][self.SPROB]) 
              if affectlen == None: #infinity time, not infected(or exposed in seir)
                 continue
              abstime = curtime + affectlen
              if smodel in ["sir","seir"] and abstime >= rtimes[minnode]:
                 continue
              if smodel == "sis":
                 minstime = [elem for elem in stimes[minnode] if elem > curtime]
                 assert len(minstime) == 1
                 if abstime >= minstime[0]:
                    continue 
              if smodel in ["si","sir","sis"]:
                 iactive.setdefault(node,abstime)   
                 if abstime < iactive[node]:
                    iactive[node] = abstime
              elif smodel == "seir":
                 eactive.setdefault(node,abstime)  
                 if abstime < eactive[node]:
                    eactive[node] = abstime
   
       if smodel in ["sir","seir","si"]:
          for times in [etimes,itimes,rtimes]:
              for node in [node for node in times.keys() if times[node] >= self.MAXSPREADTIME]:
                  del times[node]
       elif smodel == "sis":
          for node in stimes.keys():
              removetimes = [time for time in stimes[node] if time >= self.MAXSPREADTIME]
              stimes[node] -= set(removetimes)
          for node in itimes.keys():
              removetimes = [time for time in itimes[node] if time >= self.MAXSPREADTIME]
              itimes[node] -= set(removetimes)
       return stimes,etimes,itimes,rtimes
