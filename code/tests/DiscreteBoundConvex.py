import random
import unittest
import sys
import networkx as nx

class DiscreteGreedyTest(unittest.TestCase):
    """ Discrete Greedy history reconstruction  test class
    """
    
    def setUp(self,trace,timefracs,timecounts,smodel,G,noise,noisetype,infermode):
       """sets up test
         Args:
            trace: trace data
            smodel: spreading model
            distinfo: distribution information
       """
       self.G = G
       self.trace = trace  
       self.smodel = smodel
       self.timecounts = timecounts
       self.timefracs = timefracs
       self.infermode = infermode
       self.noise = noise
       self.noisetype = noisetype
       
    def testHistRecon(self):
       """ test history reconstruction
       """
       return 
        
         
    def testSiTrace(self,distinfo):
       """tests related to si traces
         Args:
            nodestate: 
       """ 
       allnodes = self.trace.keys()
       itimes = {node: self.trace[node][INFECTED] for node in allnodes if self.trace[node].has_key(INFECTED)}
       
    def testSirTrace(self,distinfo):
       """tests sir trace data
         Args:
            distinfo: distribution info
       """ 
       allnodes = self.trace.keys()
       itimes = {node: self.trace[node][INFECTED] for node in allnodes if self.trace[node].has_key(INFECTED)}
       rtimes = {node: self.trace[node][RECOVERED] for node in allnodes if self.trace[node].has_key(RECOVERED)}
       for node in rtimes.keys():
           assert rtimes[node] >= itimes[node]
       
    def testSeirTrace(self,distinfo):
       """tests seir trace data
         Args:
            distinfo: distribution info
       """ 
       allnodes = self.trace.keys()
       etimes = {node: self.trace[node][EXPOSED] for node in allnodes if self.trace[node].has_key(EXPOSED)}
       itimes = {node: self.trace[node][INFECTED] for node in allnodes if self.trace[node].has_key(INFECTED)}
       rtimes = {node: self.trace[node][RECOVERED] for node in allnodes if self.trace[node].has_key(RECOVERED)}
       for node in rtimes.keys():
           if itimes.has_key(node):
              assert rtimes[node] >= itimes[node] 
           if etimes.has_key(node):
              assert itimes[node] >= etimes[node]
              
    def testNodeStateConversion(self,trace):
       """tests node state conversion of
         Args:
            nodestate: node state
       """
       allnodes = self.trace.keys()
       itimes = {node: self.trace[node][INFECTED] for node in allnodes if self.trace[node].has_key(INFECTED)}
       stimes = {node: self.trace[node][SUSCEPTIBLE] for node in allnodes if self.trace[node].has_key(SUSCEPTIBLE)}
       alltimes = set([time for node in itimes.keys() for time in itimes[node]]).union(set([time for node in stimes.keys() for time in stimes[node]])) 
       nodestate = {}
       maxtime = sorted(alltimes)[-1]    
       nodestate = {node: {time: {SUSCEPTIBLE: -1, INFECTED: -1}}  for node in allnodes for time in alltimes}
       for node in allnodes: 
           if not stimes.has_key(node):
              assert len(itimes[node]) == 1
              nodestate[node] = {time: {SUSCEPTBLE: 0, INFECTED: 1}  for node in allnodes for time in alltimes}
              continue
           if not itimes.has_key(node):
              assert len(stimes[node]) == 1
              nodestate[node] = {time: {SUSCEPTBLE: 1, INFECTED: 0}  for node in allnodes for time in alltimes}
              continue
           nodeitimes = set(itimes[node])
           nodestimes = set(stimes[node])
           nodealltimes = sorted(list(set(nodeitimes).union(set(nodestimes))))
           for timeindex in range(1,len(nodealltimes)):
               prevtime, curtime = nodealltimes[timeindex-1:timeindex]
               for time in alltimes:
                   if time < curtime and time >= prevtime:
                      if curtime in nodeitimes: 
                         nodestate[node][time] = {INFECTED:0 , SUSCEPTIBLE = 1}
                      elif curtime in nodestimes:
                         nodestate[node][time] = {INFECTED:1 , SUSCEPTIBLE = 0}
           nodemaxtime = nodealltimes[-1]
           for time in alltimes:
               if time >= nodemaxtime and time <= maxtime:
                  if nodemaxtime in nodeitimes: 
                     nodestate[node][time] = {INFECTED:1 , SUSCEPTIBLE = 0}
                  elif nodemaxtime in nodestimes:
                     nodestate[node][time] = {INFECTED:0 , SUSCEPTIBLE = 1}

       for node in nodestate.keys():
           for time in nodestate[node].keys():
               assert nodestate[node][time][INFECTED] != -1 and nodestate[node][time][SUSCEPTIBLE] != -1

    def testSisOrder(self):
       """tests sis orderings of loopy states
       """ 
       allnodes = self.trace.keys()
       itimes = {node: self.trace[node][INFECTED] for node in allnodes if self.trace[node].has_key(INFECTED)}
       stimes = {node: self.trace[node][SUSCEPTIBLE] for node in allnodes if self.trace[node].has_key(SUSCEPTIBLE)}
       alltimes = set([time for node in itimes.keys() for time in itimes[node]]).union(set([time for node in stimes.keys() for time in stimes[node]]))
       for node in itimes.keys():
           if not stimes.has_key(node):
              assert len(itimes[node]) == 1
              continue
           nodetimes = sorted(set(itimes[node]).union(stimes[node]))
           lastflag = -1
           for time in nodetimes:
               if time in itimes[node]:
                  assert lastflag != 1 
                  lastflag = 1
               elif time in stimes[node]:
                  assert lastflag != 2 
                  lastflag = 2
        
    def testSisTrace(self,distinfo):
       """sis trace data test function
       Args:
          distinfo: dist info
       """
       #avgIcount = sum([len(itimes)[node] for node in allnodes if itimes.has_key(node)]) / float(len(itimes.keys()))
       #avgScount = sum([len(stimes)[node] for node in allnodes if stimes.has_key(node)]) / float(len(stimes.keys())) 
       testSisOrder()
       testNodeStateConversion()
       return nodestate 
    
    def testTrace(self,distinfo):
       funcname = "test{0}Trace".format(self.smodel)
       method = getAttr(self,funcname)
       method(distinfo)
    
if __name__ == '__main__':
    unittest.main()
