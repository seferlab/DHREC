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
from copy import deepcopy
from DiscreteDistribution import DiscreteDistribution as Dist
from Trace import Trace

class DiscreteTrace(Trace):
    """ Discrete Trace Generation class -> inherits from Trace
    """
    def __init__(self):
      return

    @classmethod
    def genTraces(self,G,startnodecount,evol,tracecount,smodel,modelparams={},startnode=None):
       """generates discrete traces
       Args:
          G: Graph with attributes
          startnodecount: start nodes
          evol: static/dynamic
          tracecount: number of traces
          smodel: spreading model
          modelparams: model parameters
          startnode: trace start node
       Returns:
          traces: generated traces 
       """
       if evol == "static":  
          allnodes = G.nodes()
       elif evol == "dynamic":
          self.MAXSPREADTIME = max(G.keys()) + 1
          allnodes = set([node for time in G.keys() for node in G[time].nodes()])
       traces = []
       assert smodel in ["si","sir","seir","sis"]
       for index in xrange(tracecount):
           methodname = "gen{0}Trace".format(smodel.capitalize())
           method = getattr(self,methodname)
           if startnode == None:
              random.shuffle(allnodes)
              startnodes = allnodes[0:startnodecount]
           else:
              startnodes = [int(startnode)] 
           (stimes,etimes,itimes,rtimes) = method(G,startnodes,evol,modelparams)
           trace = self.convertTimes2Trace(smodel,stimes,etimes,itimes,rtimes)
           traces.append(deepcopy(trace))
       return traces

    @classmethod
    def genSiTrace(self,G,startnodes,evol,modelparams=None):
       """generates discrete si trace 
       Args:
          G: Graph with attributes
          startnodes: start nodes
          evol: static/dynamic
       Returns:
          trace: generated trace
       """
       if evol == "static":  
          allnodes = G.nodes()
       elif evol == "dynamic":
          self.MAXSPREADTIME = max(G.keys()) + 1
          allnodes = set([node for time in G.keys() for node in G[time].nodes()]) 
       nons2ipertime = {(node1,node2): Dist.genPartDist(G[node1][node2][self.S2I][0],G[node1][node2][self.S2I][1],"reverseCdf",G[node1][node2][self.SPROB]) for node1,node2 in G.edges()}
       stimes, etimes, rtimes = {}, {}, {}
       itimes = {node : 0 for node in startnodes}
       curstate = {node: self.SUSCEPTIBLE for node in allnodes}
       for node in startnodes:
           curstate[node] = self.INFECTED
       nextstate = {node: None for node in allnodes}
       fixedcount = 0
       lastcount = len(itimes.keys()) + 1
       for time in xrange(self.MAXSPREADTIME):
           if len(itimes.keys()) == lastcount:
              fixedcount += 1
           else:
              fixedcount = 0 
           if fixedcount == 10:
              break 
           lastcount = len(itimes.keys())
           for node in allnodes:
               if curstate[node] == self.SUSCEPTIBLE:
                  if evol == "static":
                     neighset = set(spreader for spreader in G.predecessors(node) if curstate[spreader] == self.INFECTED) 
                  elif evol == "dynamic":
                     neighset = set(spreader for spreader in G[time].predecessors(node) if curstate[spreader] == self.INFECTED)
                  psus = reduce(lambda x, y: x*y, [1.0] + [nons2ipertime[(spreader,node)][time - itimes[spreader] + 1] for spreader in neighset if time - itimes[spreader] + 1 <= 10])
                  if random.random() <= (1.0 - psus):
                     nextstate[node] = self.INFECTED
                     assert node not in itimes.keys()
                     itimes[node] = time + 1
                  else:
                     nextstate[node] = self.SUSCEPTIBLE
               else:
                  nextstate[node] = self.INFECTED  
           curstate = dict(nextstate)
       return (stimes,etimes,itimes,rtimes)

    @classmethod
    def genSisTrace(self,G,startnodes,evol,modelparams):
       """generates discrete sis trace 
       Args:
          G: Graph with attributes
          startnodes: start nodes
          allnodes: all nodes
          evol: static/dynamic
       Returns:
          trace: generated trace
       """
       self.MAXSPREADTIME = modelparams["sismaxtime"] 
       if evol == "static":  
          allnodes = G.nodes()
       elif evol == "dynamic":
          self.MAXSPREADTIME = max(G.keys()) + 1
          allnodes = set([node for time in G.keys() for node in G[time].nodes()])   
       nons2ipertime = {(node1,node2): Dist.genPartDist(G[node1][node2][self.S2I][0],G[node1][node2][self.S2I][1],"reverseCdf",G[node1][node2][self.SPROB]) for node1,node2 in G.edges()}
       noni2spertime = {node: Dist.genPartDist(G.node[node][self.I2S][0],G.node[node][self.I2S][1],"reverseCdf") for node in G.nodes()}
       itimes = {node : [] for node in allnodes}
       for node in startnodes:
           itimes[node] = [0]
       stimes = {node : [0] for node in allnodes}
       for node in startnodes:
           stimes[node] = []
       etimes, rtimes = {}, {}
       curstate = {node: self.SUSCEPTIBLE for node in allnodes}
       for node in startnodes:
           curstate[node] = self.INFECTED
       nextstate = {node: None for node in allnodes}
       inodes = set(startnodes)
       for time in xrange(self.MAXSPREADTIME):
           if len(inodes) == 0:
              break
           for node in allnodes:
               if curstate[node] == self.SUSCEPTIBLE:
                  if evol == "static":
                     neighset = set(spreader for spreader in G.predecessors(node) if curstate[spreader] == self.INFECTED) 
                  elif evol == "dynamic":
                     neighset = set(spreader for spreader in G[time].predecessors(node) if curstate[spreader] == self.INFECTED)
                  psus = reduce(lambda x, y: x*y, [1.0] + [nons2ipertime[(spreader,node)][time - itimes[spreader][-1] + 1] for spreader in neighset if time - itimes[spreader][-1] + 1 <= 10])
                  if random.random() <= (1.0 - psus):
                     assert node not in inodes
                     inodes.add(node)
                     nextstate[node] = self.INFECTED
                     itimes[node].append(time + 1)
                  else:
                     nextstate[node] = self.SUSCEPTIBLE
               elif curstate[node] == self.INFECTED:
                  pinfect = 0.0
                  timedif = time - itimes[node][-1] + 1 
                  if timedif <= 10:
                     pinfect = noni2spertime[node][timedif]
                  if random.random() <= (1.0 - pinfect):
                     assert node in inodes
                     inodes.remove(node)
                     nextstate[node] = self.SUSCEPTIBLE
                     stimes[node].append(time + 1)
                  else:
                     nextstate[node] = self.INFECTED
           curstate = dict(nextstate)    
       return (stimes,etimes,itimes,rtimes)
   
    @classmethod
    def genSirTrace(self,G,startnodes,evol,modelparams=None):
       """generates discrete sir trace 
       Args:
          G: Graph with attributes
          startnodes: start nodes
          evol: static/dynamic
       Returns:
          trace: generated trace
       """
       if evol == "static":  
          allnodes = G.nodes()
       elif evol == "dynamic":
          self.MAXSPREADTIME = max(G.keys()) + 1
          allnodes = set([node for time in G.keys() for node in G[time].nodes()])  
       nons2ipertime = {(node1,node2): Dist.genPartDist(G[node1][node2][self.S2I][0],G[node1][node2][self.S2I][1],"reverseCdf",G[node1][node2][self.SPROB]) for node1,node2 in G.edges()}
       noni2rpertime = {node: Dist.genPartDist(G.node[node][self.I2R][0],G.node[node][self.I2R][1],"reverseCdf") for node in G.nodes()}
       stimes, etimes = {}, {}
       itimes = {node : 0 for node in startnodes}
       rtimes = {}
       curstate = {node: self.SUSCEPTIBLE for node in allnodes}
       for node in startnodes:
           curstate[node] = self.INFECTED
       nextstate = {node: None for node in allnodes}
       for time in xrange(self.MAXSPREADTIME):
           flag = False
           for node in allnodes:
               if curstate[node] == self.INFECTED:
                  flag = True
                  break
           if not flag:
              break
           for node in allnodes:
               if curstate[node] == self.SUSCEPTIBLE:
                  if evol == "static":
                     neighset = set(spreader for spreader in G.predecessors(node) if curstate[spreader] == self.INFECTED) 
                  elif evol == "dynamic":
                     neighset = set(spreader for spreader in G[time].predecessors(node) if curstate[spreader] == self.INFECTED)
                  psus = reduce(lambda x, y: x*y, [1.0] + [nons2ipertime[(spreader,node)][time - itimes[spreader] + 1] for spreader in neighset if time - itimes[spreader] + 1 <= 10])
                  if random.random() <= (1.0 - psus):
                     nextstate[node] = self.INFECTED
                     assert node not in itimes.keys()
                     itimes[node] = time + 1
                  else:
                     nextstate[node] = self.SUSCEPTIBLE
               elif curstate[node] == self.INFECTED:
                  pinfect = 0.0
                  timedif = time - itimes[node] + 1
                  if timedif <= 10:
                     pinfect = noni2rpertime[node][timedif]
                  if random.random() <= (1.0 - pinfect):
                     nextstate[node] = self.RECOVERED
                     assert node not in rtimes.keys() and node in itimes.keys()
                     rtimes[node] = time + 1
                  else:
                     nextstate[node] = self.INFECTED
               else:   
                  nextstate[node] = self.RECOVERED
           curstate = dict(nextstate)
       return (stimes,etimes,itimes,rtimes)
   
    @classmethod
    def genSeirTrace(self,G,startnodes,evol,modelparams=None):
       """generates discrete seir trace 
       Args:
          G: Graph with attributes
          startnodes: start nodes
          evol: static/dynamic
       Returns:
          trace: generated trace
       """
       if evol == "static":  
          allnodes = G.nodes()
       elif evol == "dynamic":
          self.MAXSPREADTIME = max(G.keys()) + 1
          allnodes = set([node for time in G.keys() for node in G[time].nodes()])  
       nons2epertime = {(node1,node2): Dist.genPartDist(G[node1][node2][self.S2E][0],G[node1][node2][self.S2E][1],"reverseCdf",G[node1][node2][self.SPROB]) for node1,node2 in G.edges()}
       noni2rpertime = {node: Dist.genPartDist(G.node[node][self.I2R][0],G.node[node][self.I2R][1],"reverseCdf") for node in G.nodes()}
       none2ipertime = {node: Dist.genPartDist(G.node[node][self.E2I][0],G.node[node][self.E2I][1],"reverseCdf") for node in G.nodes()}
       itimes = {node : 0 for node in startnodes}
       stimes, rtimes,etimes = {}, {}, {}
       curstate = {node: self.SUSCEPTIBLE for node in allnodes}
       for node in startnodes:
           curstate[node] = self.INFECTED
       nextstate = {node: None for node in allnodes}
       for time in xrange(self.MAXSPREADTIME):
           flag = False
           for node in allnodes:
               if curstate[node] in [self.INFECTED, self.EXPOSED]:
                  flag = True
                  break
           if not flag:
              break 
           for node in allnodes:
               if curstate[node] == self.SUSCEPTIBLE:
                  if evol == "static":
                     neighset = set(spreader for spreader in G.predecessors(node) if curstate[spreader] == self.INFECTED) 
                  elif evol == "dynamic":
                     neighset = set(spreader for spreader in G[time].predecessors(node) if curstate[spreader] == self.INFECTED)
                  psus = reduce(lambda x, y: x*y, [1.0] + [nons2epertime[(spreader,node)][time - itimes[spreader] + 1] for spreader in neighset if time - itimes[spreader] + 1 <= 10])
                  if random.random() <= (1.0 - psus):
                     nextstate[node] = self.EXPOSED
                     assert node not in etimes.keys()
                     etimes[node] = time + 1
                  else:
                     nextstate[node] = self.SUSCEPTIBLE
               elif curstate[node] == self.EXPOSED:
                  pexposed = 0.0
                  timedif = time - etimes[node] + 1
                  if timedif <= 10: 
                     pexposed = none2ipertime[node][timedif]
                  if random.random() <= (1.0 - pexposed):
                     nextstate[node] = self.INFECTED
                     assert node not in itimes.keys() and node in etimes.keys()
                     itimes[node] = time + 1
                  else:
                     nextstate[node] = self.EXPOSED  
               elif curstate[node] == self.INFECTED:
                  pinfect = 0.0
                  timedif = time - itimes[node] + 1
                  if timedif <= 10: 
                     pinfect = noni2rpertime[node][timedif]
                  if random.random() <= (1.0 - pinfect):
                     nextstate[node] = self.RECOVERED
                     assert node not in rtimes.keys() and node in itimes.keys()
                     rtimes[node] = time + 1
                  else:
                     nextstate[node] = self.INFECTED
               else:   
                  nextstate[node] = self.RECOVERED
           curstate = dict(nextstate)
       return (stimes, etimes, itimes, rtimes)      




