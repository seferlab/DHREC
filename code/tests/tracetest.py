import random
import unittest
import sys
sys.path.append("../lib")
sys.path.append("..")
import networkx as nx
import myutilities as myutil
import os
import math
import itertools
import operator
import autohistorytracegen
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace
from DiscreteDistribution import DiscreteDistribution as DisDist
from ContinuousDistribution import ContinuousDistribution as ContDist
from Score import Score

class TraceTest(unittest.TestCase):
    """ Trace data test class
    """
    def __init__(self,testname,traces=None,smodel=None,distinfo=None,G=None,noise=None,samplerate=None,prob=None):
       """constructor
         Args:
            testname: test to run
            traces: trace datas dictionary
            smodel: spreading model
            distinfo: distribution information
            G: graph
            noise:
            samplerate:
            prob:
       """
       super(TraceTest, self).__init__(testname)
       self.traces = traces
       self.smodel = smodel
       self.distinfo = distinfo
       self.noise = noise
       self.samplerate = samplerate
       self.prob = prob
       self.G = G
       self.isTraceEmpty()

    def testTraceValidity(self):
        """ tests whether trace is valid
        """
        for trace in self.traces:
            if not Trace.IsTraceNoisy([trace]):
               for node in trace.keys():
                   alltimes = trace[node].values()
                   self.assertEqual(len(alltimes),len(set(alltimes)))
            else:
               alltimes = []
               for node in trace.keys():
                   alltimes.append(trace[node].keys())
                   for time in trace[node].keys():
                       self.assertLess(abs(sum(trace[node][time].values())- 1.0), 0.01)
               for index1 in xrange(1,len(alltimes)):
                   self.assertEqual(len(set(alltimes[index1-1]).difference(set(alltimes[index1]))) , 0)
                   self.assertEqual(len(set(alltimes[index1]).difference(set(alltimes[index1-1]))) , 0)

    def testNoisyTrace(self):
        """ tests:
        1- they add up to 1
        2- noisedegree and current state makes sense
        3- if allnodes have the same time steps
        """
        self.testTraceValidity() 
        for trace in self.traces:
            self.assertTrue(Trace.IsTraceNoisy([trace]))
            difsum,count = 0.0,0
            for node in trace.keys():
                for time in trace[node].keys():
                    for state in trace[node][time].keys():
                        difsum += min(1.0 - trace[node][time][state],trace[node][time][state])
                        count += 1
            self.assertLess(abs(self.noise - difsum/float(count)), 1.0)      
                          
    def testSampledTrace(self):
        """tests subsampling method
           1- they must obey samplerate
        """
        self.testTraceValidity() 
        for trace in self.traces:
            alltimes = Trace.getAllTimes([trace])
            for time in alltimes:
                assert (time % self.samplerate) == 0

    def testCoverValidity(self,trace):
        """tests validity of timing between states
        Args:
           trace:
        Returns:
           bool: 
        """  
        allnodes = trace.keys()  
        if self.smodel == "si":
           sentimes = {node: trace[node][Trace.INFECTED] for node in allnodes if trace[node].has_key(Trace.INFECTED)}
           rectimes = sentimes
        elif self.smodel == "sir":
           sentimes = {node: trace[node][Trace.INFECTED] for node in allnodes if trace[node].has_key(Trace.INFECTED)}
           rectimes = sentimes
           rtimes = {node: trace[node][Trace.RECOVERED] for node in allnodes if trace[node].has_key(Trace.RECOVERED)}
        elif self.smodel == "seir":    
           sentimes = {node: trace[node][Trace.INFECTED] for node in allnodes if trace[node].has_key(Trace.INFECTED)}
           rectimes = {node: trace[node][Trace.EXPOSED] for node in allnodes if trace[node].has_key(Trace.EXPOSED)}
           rtimes = {node: trace[node][Trace.RECOVERED] for node in allnodes if trace[node].has_key(Trace.RECOVERED)}
        senpairs = sorted(sentimes.iteritems(), key=operator.itemgetter(1))
        recpairs = sorted(rectimes.iteritems(), key=operator.itemgetter(1))
        startnodes = []
        for node,time in senpairs:
            if time != 0: 
               break 
            startnodes.append(node)
        avgdur,trans = 0.0, 0
        for index1 in xrange(len(recpairs)):
            node1, itime1 = recpairs[index1]
            if node1 in startnodes:
               continue 
            edgeexist = False
            for index2 in xrange(len(senpairs)):
                node2, itime2 = senpairs[index2]
                if itime2 < itime1 and self.G.has_edge(node2,node1): 
                   if self.smodel in ["sir","seir"] and rtimes.has_key(node2) and rtimes[node2] < itime1:
                      continue 
                   edgeexist = True
                   avgdur += (itime1 - itime2)
                   trans += 1 
                   break
            self.assertTrue(edgeexist)
        if trans != 0:
           avgdur /= float(trans) 
        return True
               
    def testSiTrace(self):
       """si trace generation test
          This test verifies:
          1- startnode existence
          2- Edge existence between consecutive time points
          3- Compatibility of infection times and average infected nodes with given parameters(avgi,avgdur)
       """
       avgi = 0.0
       for trace in self.traces:
           allnodes = trace.keys()
           itimes = {node: trace[node][Trace.INFECTED] for node in allnodes if trace[node].has_key(Trace.INFECTED)}
           avgi += len(itimes.keys())
           self.assertTrue(self.testCoverValidity(trace))
       avgi /= float(len(self.traces))

    def testSirTrace(self):
        """sir trace generation test
           This test verifies:
           1- startnode existence
           2- Edge existence between consecutive time points
           3- Compatibility of infection times and recovery times with given parameters
           4- recovery time must be greater than infection time
        """
        avgi, avgr = 0.0, 0.0 
        for trace in self.traces:   
            allnodes = trace.keys()
            itimes = {node: trace[node][Trace.INFECTED] for node in allnodes if trace[node].has_key(Trace.INFECTED)}
            rtimes = {node: trace[node][Trace.RECOVERED] for node in allnodes if trace[node].has_key(Trace.RECOVERED)}
            self.assertTrue(self.testCoverValidity(trace))
            avgrtime,rcount = 0.0, 0
            for node in allnodes:
                if itimes.has_key(node) and rtimes.has_key(node):
                   self.checkGreater(rtimes[node],itimes[node])
                   rcount += 1
                   avgrtime += rtimes[node] - itimes[node]
            avgrtime /= float(rcount)

            if self.prob == "cont":
               if rcount < 10 and self.samplerate != 0:
                  continue
               realavg = ContDist.getDistMean(self.distinfo[Trace.I2R][0],self.distinfo[Trace.I2R][1])
               self.assertLess(abs(avgrtime - realavg) , realavg / 2.0)
            elif self.prob == "dis":
               pass 

    def checkGreater(self,val1,val2):
        """check first time is before second time
           equality is not accepted when samplerate =0
        """
        if self.samplerate == 0:
           self.assertGreater(val1,val2)
        else:     
           self.assertGreaterEqual(val1,val2)

    def isTraceEmpty(self):
        """ checks whether trace is empty(not valid)
        """
        for trace in self.traces:
            if Trace.IsTraceNoisy([trace]):
               if len(trace.keys()) == 0:
                  return False
               for node in trace.keys():
                   if len(trace[node].keys()) == 0:
                      return False
                   for time in trace[node].keys():
                       if len(trace[node][time].keys()) == 0:
                          return False
               return True        
            else:
               if len(trace.keys()) == 0:
                  return False
                  for node in trace.keys():
                      if len(trace[node].keys()) == 0: 
                          return False
                  return True
               
    def testSeirTrace(self):
        """seir trace generation test
           This test verifies:
           1- startnodes from filename
           2- Edge existence between consecutive time points
           3- Compatibility of infection,exposed,recovery times with given parameters
           4- recovery time >= infection time and infection time >= exposed time 
        """
        for trace in self.traces:
            allnodes = trace.keys()
            etimes = {node: trace[node][Trace.EXPOSED] for node in allnodes if trace[node].has_key(Trace.EXPOSED)}
            itimes = {node: trace[node][Trace.INFECTED] for node in allnodes if trace[node].has_key(Trace.INFECTED)}
            rtimes = {node: trace[node][Trace.RECOVERED] for node in allnodes if trace[node].has_key(Trace.RECOVERED)}
            self.assertTrue(self.testCoverValidity(trace))
            avgrtime,rcount = 0.0, 0
            for node in allnodes:
                if itimes.has_key(node) and rtimes.has_key(node):
                   self.checkGreater(rtimes[node],itimes[node])
                   rcount += 1
                   avgrtime += rtimes[node] - itimes[node]
            avgrtime /= float(rcount)
            avgitime,icount = 0.0, 0
            for node in allnodes:
                if itimes.has_key(node) and etimes.has_key(node):
                   self.checkGreater(itimes[node],etimes[node])
                   icount += 1
                   avgrtime += itimes[node] - etimes[node]
            avgitime /= float(icount)
            
            if self.prob == "cont":
               if rcount < 20 or self.samplerate != 0 :
                  continue
               else:
                  realavg = ContDist.getDistMean(self.distinfo[Trace.I2R][0],self.distinfo[Trace.I2R][1])
                  self.assertLess(abs(avgrtime - realavg) , realavg / 2.0)
               if icount < 20 or self.samplerate != 0:
                  continue   
               else:
                  realavg = ContDist.getDistMean(self.distinfo[Trace.E2I][0],self.distinfo[Trace.E2I][1])
                  self.assertLess(abs(avgitime - realavg) , realavg / 2.0)
            elif self.prob == "dis":
               pass
            
    def testProbGen(self):
       """ tests probability generating methods
       This test verifies:
           1- Probability generation for both discrete and continous cases
       """
       if self.prob == "cont":
          sprob = self.distinfo[Trace.SPROB]
          for key in self.distinfo.keys():
              if key == Trace.SPROB:
                 continue
              self.assertLess(abs(ContDist.getContCdf(self.distinfo[key][0],self.distinfo[key][1],100000,sprob) - sprob), 0.01)

              #random number gen test
              numcount= 5000
              avgnum = sum([ContDist.genContRandNum(self.distinfo[key][0],self.distinfo[key][1]) for index in xrange(numcount)]) / float(numcount)
              self.assertLess(abs(avgnum - ContDist.getDistMean(self.distinfo[key][0],self.distinfo[key][1])) , avgnum/4)  
              intlen = 0.01
              xlist = range(1,7)
              combs = [] 
              for index1 in xrange(len(xlist)):
                  for index2 in xrange(index1+1, len(xlist)):
                      combs.append((xlist[index1],xlist[index2]))
              for lower,upper in combs:
                  area = 0.0
                  for index in xrange(1,int((upper-lower)/intlen)):
                      x1 = lower + ((index-1) * intlen)
                      x2 = lower + (index * intlen)
                      y1 = ContDist.getContPdf(self.distinfo[key][0],self.distinfo[key][1], x1, sprob)
                      y2 = ContDist.getContPdf(self.distinfo[key][0],self.distinfo[key][1], x2, sprob)
                      area += intlen * ((y1 + y2) / 2.0)
                  realarea = ContDist.getContCdf(self.distinfo[key][0],self.distinfo[key][1],upper,sprob) - ContDist.getContCdf(self.distinfo[key][0],self.distinfo[key][1],lower,sprob)
                  self.assertLess(area,1.0)
                  self.assertLess(realarea,1.0)
                  self.assertLess(abs(area - realarea) , 0.01)
       elif self.prob == "dis":
          sprob = self.distinfo[Trace.SPROB]
          for key in self.distinfo.keys():
              if key != Trace.SPROB:
                 methodname = "get{0}Pdf".format(self.distinfo[key][0].capitalize())
                 method = getattr(DisDist,methodname)
                 self.assertLess(abs(1.0 - sum(method(self.distinfo[key][1],10).values())), 0.01)
                 ratios = DisDist.genPartDist(self.distinfo[key][0],self.distinfo[key][1],"normal",sprob,20)
                 probs,restmul = {}, 1.0
                 for time in sorted(ratios.keys()):
                     probs[time] = ratios[time] * restmul
                     restmul *= ratios[time]
                 self.assertLess(abs(sprob - sum(probs.values())), 0.01)
                 ratios = DisDist.genPartDist(self.distinfo[key][0],self.distinfo[key][1],"normalCdf",sprob,20)
                 restmul = reduce(operator.mul,ratios.values(),1.0)
                 self.assertLess(abs(sprob - restmul), 0.02)
                 ratios = DisDist.genPartDist(self.distinfo[key][0],self.distinfo[key][1],"reverseCdf",sprob,20)
                 restmul = reduce(operator.mul,ratios.values(),1.0)
                 self.assertLess(abs(1.0 - sprob - restmul), 0.02)
              
    def testSisStateConversion(self):
        """ tests sis node state conversion
            This test verifies:
            1- We don't lose info by changing format during noise addition 
            2- Noise addition and subsampling results are valid
        """
        for trace in self.traces:
            allnodes = trace.keys()
            itimes = {node: trace[node][Trace.INFECTED] for node in allnodes if trace[node].has_key(Trace.INFECTED)}
            stimes = {node: trace[node][Trace.SUSCEPTIBLE] for node in allnodes if trace[node].has_key(Trace.SUSCEPTIBLE)}
            alltimes = set([time for node in itimes.keys() for time in itimes[node]]).union(set([time for node in stimes.keys() for time in stimes[node]]))
            nodestate = {}
            maxtime = sorted(alltimes)[-1]    
            for node in allnodes:
                if not stimes.has_key(node):
                   self.assertEqual(len(itimes[node]),1)
                   nodestate[node] = {time: {Trace.SUSCEPTIBLE: 0, Trace.INFECTED: 1}  for node in allnodes for time in alltimes}
                   continue
                if not itimes.has_key(node):
                   self.assertEqual(len(stimes[node]),1)
                   nodestate[node] = {time: {Trace.SUSCEPTIBLE: 1, Trace.INFECTED: 0}  for node in allnodes for time in alltimes}
                   continue
                nodetimes = sorted(list(set(stimes[node]).union(set(itimes[node]))))
                for timeindex in xrange(1,len(nodetimes)):
                    prevtime, curtime = nodetimes[timeindex-1 : timeindex]
                    for time in alltimes:
                        if time < curtime and time >= prevtime:
                           if curtime in itimes[node]: 
                              nodestate[node][time] = {Trace.INFECTED:0 , Trace.SUSCEPTIBLE: 1}
                           elif curtime in stimes[node]:
                              nodestate[node][time] = {Trace.INFECTED:1 , Trace.SUSCEPTIBLE: 0}
                maxnodetime = nodetimes[-1]
                for time in alltimes:
                    if time >= maxnodetime and time <= maxtime:
                       if maxnodetime in itimes[node]: 
                          nodestate[node][time] = {Trace.INFECTED:1 , Trace.SUSCEPTIBLE: 0}
                       elif maxnodetime in stimes[node]:
                          nodestate[node][time] = {Trace.INFECTED:0 , Trace.SUSCEPTIBLE: 1}
        for node in nodestate.keys():
            for time in nodestate[node].keys():
                self.assertEqual(nodestate[node][time][Trace.INFECTED] + nodestate[node][time][Trace.SUSCEPTIBLE], 1)
                
    def testSisOrder(self):
       """tests time ordering of nodes state transition times in sis model
       """
       for trace in self.traces():   
           allnodes = trace.keys()
           itimes = {node: trace[node][Trace.INFECTED] for node in allnodes if trace[node].has_key(Trace.INFECTED)}
           stimes = {node: trace[node][Trace.SUSCEPTIBLE] for node in allnodes if trace[node].has_key(Trace.SUSCEPTIBLE)}
           alltimes = set([time for node in itimes.keys() for time in itimes[node]]).union(set([time for node in stimes.keys() for time in stimes[node]]))
           for node in itimes.keys():
               if not stimes.has_key(node):
                  self.assertEqual(len(itimes[node]), 1)
                  continue
               nodetimes = sorted(set(itimes[node]).union(stimes[node]))
               lastflag = -1
               for time in nodetimes:
                   if time in itimes[node]:
                      self.assertNotEqual(lastflag, 1) 
                      lastflag = 1
                   elif time in stimes[node]:
                      self.assertNotEqual(lastflag, 2)  
                      lastflag = 2
        
    def testSisTrace(self):
       """ test sis generated trace data
       """
       testSisOrder()
       testSisStateConversion()


def getManyDist(smodel,prob):
    if prob == "cont": 
       dists = [("expo",(1.0,)),("expo", (2.0,)),("expo", (3.0,)),("expo", (4.0,)),("expo", (5.0,)),("expo",(0.5,)),("weibull", (2.5, 9.5)) ,("weibull", (1.5, 4.5)) ,("weibull", (2.5, 4.5)),("weibull", (3.5, 5.5)), ("rayleigh",(0.5,)),("rayleigh",(1.0,)),("rayleigh",(1.5,)),("rayleigh",(2.0,))]
    elif prob == "dis":
       dists = [("expo",(1.0,)),("expo", (2.0,)),("expo", (3.0,)),("expo", (4.0,)),("expo", (5.0,)),("expo",(0.5,)),("weibull", (2.5, 9.5)) ,("weibull", (1.5, 4.5)) ,("weibull", (2.5, 4.5)),("weibull", (3.5, 5.5)), ("rayleigh",(0.5,)),("rayleigh",(1.0,)),("rayleigh",(1.5,)),("rayleigh",(2.0,))]
    sprobs = [0.05,0.08,0.1,0.2,0.25,0.5,0.75]
    random.shuffle(dists)
    random.shuffle(sprobs)
    spreadparams = []
    if smodel == "si":
       for s2i,sprob in itertools.product(dists,sprobs):
           sparam = {Trace.S2I: s2i, Trace.SPROB: sprob}
           spreadparams.append(sparam)
    elif smodel == "sir":
       for s2i,i2r,sprob in itertools.product(dists,dists,sprobs):
           sparam = {Trace.S2I: s2i, Trace.I2R: i2r, Trace.SPROB: sprob}
           spreadparams.append(sparam)
    elif smodel == "seir":
       for s2e,e2i,i2r,sprob in itertools.product(dists,dists,dists,sprobs):
           sparam = {Trace.S2E: s2e, Trace.E2I: e2i, Trace.I2R: i2r, Trace.SPROB: sprob}
           spreadparams.append(sparam)
    elif smodel == "sis":
       for s2i,i2s,sprob in itertools.product(dists,dists,sprobs):
           sparam = {Trace.S2I: s2i, Trace.I2S: i2s,Trace.SPROB: sprob}
           spreadparams.append(sparam)
    return spreadparams

  
def getTests(smodel,noise,samplerate,prob,G,distinfo,traces):
    """ generates tests
    """
    suiteFew = unittest.TestSuite()
    suiteFew.addTest(TraceTest("testProbGen",traces,smodel,distinfo,G,noise,samplerate,prob))
    
    #if noise != 0.0:
    #   suiteFew.addTest(TraceTest("testNoisyTrace",traces,smodel,distinfo,G,noise,samplerate,prob))
    #   return suiteFew
    
    if samplerate != 0:
       suiteFew.addTest(TraceTest("testSampledTrace",traces,smodel,distinfo,G,noise,samplerate,prob))
    if smodel == "si":
       suiteFew.addTest(TraceTest("testSiTrace",traces,smodel,distinfo,G,noise,samplerate,prob))
    elif smodel == "sir":   
       suiteFew.addTest(TraceTest("testSirTrace",traces,smodel,distinfo,G,noise,samplerate,prob))
    elif smodel == "sis":   
       suiteFew.addTest(TraceTest("testSisTrace",traces,smodel,distinfo,G,noise,samplerate,prob))
    elif smodel == "seir":
       suiteFew.addTest(TraceTest("testSeirTrace",traces,smodel,distinfo,G,noise,samplerate,prob))
    return suiteFew

def main():
    """ main function for testing traces
    """
    noisetype = "StateChange"
    samplerates = [0]
    noises = [0.0,0.1,0.2,0.3,0.4,0.5,0.7,0.75,0.9,0.95]
    realdatas = ["real","syn"]
    evols = ["static"]
    #tracecounts = [5,10,15,20,25,30,50,60,75,100,150,200]
    tracecounts = [5]
    startcounts = [1,2,3,4,5,6,7,8,9,10,11,12]
    probs = ["dis"]
    smodels = ["si","sir","seir"]
    graphpref = "graphs"
    tracepref = "testtraces"
    configpref = "testconfig"
    tracefolderpref = "."
    graphfolderpref = ".."
    configfolderpref = "."
    pbsfolder = "pbsfolder"
    filesamplecount = 10
    sismaxtime = 100
    
    parts = list(itertools.product(samplerates,noises,realdatas,evols,tracecounts,probs,smodels,startcounts))
    random.shuffle(parts)
    for (samplerate,noise,realdata,evol,tracecount,prob,smodel,startcount) in parts:
        if samplerate == 0 and noise != 0.0:
           continue  
        for spreadparam in getManyDist(smodel,prob):
            sparamstr = Trace.getSpreadFolder(smodel,spreadparam,prob)
            graphfolder,tracefolder,configfolder = autohistorytracegen.genMainFolders(graphpref,graphfolderpref,tracepref,tracefolderpref,configpref,configfolderpref,pbsfolder,realdata,smodel,evol,prob)
            path2info = autohistorytracegen.returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount,prob,smodel)
            for path in path2info.keys():
                Gname,filetracefolder = path2info[path]
                G = InputOutput.readGraphAndParams(path)
                extensionstr = "-".join(path.replace("../","").split("/")[1:])
                vararr = {}
                vararr["ongraph"] = True
                vararr["traceoutfolder"] = filetracefolder
                vararr["graphfilename"] = path
                vararr["tracecount"] = tracecount
                vararr["samplerate"] = samplerate
                vararr["startnodecount"] = startcount
                vararr["noise"] = noise
                vararr["noisetype"] = noisetype
                vararr["smodel"] = smodel
                vararr["evol"] = evol
                vararr["dists"] = spreadparam
                vararr["dist"] = prob
                if smodel == "sis":
                   vararr["sismaxtime"] = sismaxtime
                configpath = "{0}/{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}-{9}.config".format(configfolder,extensionstr,realdata,evol,smodel,prob,tracecount,sparamstr,samplerate,noise)
                autohistorytracegen.genTraceConfigFile(configpath,vararr)
                code = "python ../tracegen.py {0}".format(configpath)
                os.system(code)
                #pbsfilename = "{0}-{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}.pbs".format(realdata,evol,smodel,prob,extensionstr,tracecount,sparamstr,samplerate,noise)
                #pbsrunner(code,pbsfolder,pbsfilename)
                noisesamplestr = "{0}-{1}-{2}".format(noisetype,noise,samplerate)
                traceoutfolder = "{0}/{1}/{2}".format(filetracefolder,noisesamplestr,startcount)
                (traces,allnodes) = InputOutput.readTraces(traceoutfolder,tracecount,smodel)
                [os.system("rm -rf {0}".format(folder)) for folder in [traceoutfolder,configpath]]
                suiteFew = getTests(smodel,noise,samplerate,prob,G,spreadparam,traces)
                unittest.TextTestRunner(verbosity = 2).run(suiteFew)
                
    
if __name__ == '__main__':
    main()
