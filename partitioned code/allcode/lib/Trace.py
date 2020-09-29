import networkx as nx
import numpy as np
import scipy as sp
import random
import math
import sys
import os
import myutilities as myutil
import gzip
import cPickle
import string 
import itertools
from copy import deepcopy

class Trace():     
    """ Methods related to trace operations
    """
    def __init__():
       return

    MAXSPREADTIME = 1000000000000000
    INFECTED = "i"
    EXPOSED = "e"
    SUSCEPTIBLE = "s"
    RECOVERED = "r"
    I2R = "i2r"
    E2I = "e2i"
    S2I = "s2i"
    I2S = "i2s"
    S2E = "s2e"
    SPROB = "sprob"

    @classmethod
    def IsTraceNoisy(self,traces):
        """checks whether all traces are perfect
        Args:
          trace: trace
        Returns:
          flag: True/False(whether trace is noisy(True))
        """
        for trace in traces:
            if type(trace[trace.keys()[0]].keys()[0]) != str:
               return True 
        return False

    @classmethod
    def getAllNodes(self,traces):
        """returns all seen nodes
        Args:
          traces: all seen traces
        Returns:
          allnodes: allnodes seen in traces  
        """
        return set([node for trace in traces for node in trace.keys()])
        
    @classmethod
    def getAllTimes(self,traces):
        """returns all seen times
        Args:
          traces: all seen traces
        Returns:
          alltimes: alltimes seen in traces  
        """
        if self.IsTraceNoisy(traces):
           return set([time for trace in traces for node in trace.keys() for time in trace[node].keys()])
        else:    
           return set(itertools.chain(*[trace[node].values() for trace in traces for node in trace.keys() if len(trace[node].keys()) != 0]))

    @classmethod
    def getLseProb(self,trace,sortedtimes,smodel,inftime):
       """returns least square error probabilities
       Args:
          trace: tracedata 
          sortedtimes: alltimes sorted
          smodel: spreading model
          inftime: infinitytime needed to model uninfected probability of nodes
       Returns:
          returns the LSE infection probabilities for si model
       """
       assert smodel == "si"
       iprobs = {}
       for node in trace.keys():
           iprobs[node] = {}
           posinfected = [time for time in sortedtimes if trace[node][time][Trace.INFECTED] > 0.01] #possible infection times for node
           if len(posinfected) == 0:
              iprobs[node][inftime] = 1.0
              continue
           iprobs[node] = {time: 0.0 for time in sortedtimes if time not in posinfected}
           for itime in posinfected:
               err = sum([trace[node][time][Trace.INFECTED]**2 for time in sortedtimes if time < itime])
               err += sum([(1.0-trace[node][time][Trace.INFECTED])**2 for time in sortedtimes if time >= itime])
               if err == 0.0:
                  iprobs[node] = {} 
                  iprobs[node][itime] = 1.0
                  break 
               iprobs[node][itime] = 1.0 / err
           if err == 0.0:
              continue 
           err = sum([trace[node][time][Trace.INFECTED]**2 for time in sortedtimes])   
           iprobs[node][inftime] = 1.0 / err
           totprob = sum(iprobs[node].values())
           iprobs[node] = {key: iprobs[node][key] / float(totprob) for key in iprobs[node].keys()}
       return iprobs

    @classmethod
    def roundTraceLse(self,traces,rparam):
        """ rounds traces by lse (under categorical rounding)
        Args:
           traces: noisy traces
           rparam: params
        Returns:
           newtraces: new rounded traces
        """
        INFTIME = 10000000
        newtraces = []
        for trace in traces:
            alltimes = sorted(list(set([time for node in trace.keys() for time in trace[node].keys()])))
            node2probs = self.getLseProb(trace,alltimes,"si",INFTIME)
            newtrace = {}
            for node in trace.keys():
                probint,cursum = {}, 0.0
                for time in node2probs[node].keys():
                    probint[(cursum,cursum+node2probs[node][time])] = time
                    cursum += node2probs[node][time]
                p = random.random()
                for (left,right) in probint.keys():
                    if left <= p and p <= right:
                       roundtime = probint[(left,right)]
                       break
                if roundtime != INFTIME:    
                   newtrace[node] = {Trace.INFECTED: roundtime}
            newtraces.append(newtrace)       
        return newtraces

    @classmethod
    def roundTraceRandom(self,traces,rparam):
        """ rounds traces randomly
        Args:
           traces: noisy traces
           rparam: params
        Returns:
           newtraces: new rounded traces
        """
        newtraces = []
        for trace in traces:
            alltimes = sorted(list(set([time for node in trace.keys() for time in trace[node].keys()])))
            newtrace = {}
            for node in trace.keys():
                for time in alltimes:
                    if trace[node][time].has_key(Trace.INFECTED):
                       break
                else:
                   continue
                isum = sum([trace[node][time][Trace.INFECTED] for time in alltimes])
                normdict = {time: trace[node][time][Trace.INFECTED]/float(isum) for time in alltimes}
                probint,cursum = {}, 0.0
                for time in alltimes:
                    probint[(cursum,cursum + normdict[time])] = time
                    cursum += normdict[time]
                p = random.random()
                for (left,right) in probint.keys():
                    if left <= p and p <= right:
                       roundtime = probint[(left,right)]
                       break
                newtrace[node] = {Trace.INFECTED: roundtime}
            newtraces.append(newtrace)       
        return newtraces
    
    @classmethod  
    def roundTrace(self,traces,rmethod,rparam):
        """rounds probabilstic traces for other algos
           assumes categorical rounding as baseline
        """
        assert rmethod in ["lse","random"]
        methodname = "roundTrace{0}".format(rmethod.capitalize())
        method = getattr(self,methodname)
        return method(traces,rparam)

    @classmethod 
    def convertTimes2Trace(self,smodel,stimes,etimes,itimes,rtimes):
       """converts time to trace format for pkl
       Args:
          smodel: spreading model
          itimes, etimes, rtimes, stimes: state times 
       Returns:
          trace: trace for pkl format
       """
       trace = {}
       statemap = {0: self.SUSCEPTIBLE, 1: self.EXPOSED, 2: self.INFECTED, 3: self.RECOVERED}
       times = [stimes, etimes, itimes, rtimes]
       for index in xrange(4):
           for node in times[index].keys(): 
               trace.setdefault(node,{})
               trace[node][statemap[index]] = times[index][node]
       return trace

    @classmethod
    def parseTransitionWithDist(self,info,name):
        """ parse given transition string to distribution with dist name
        Args:
           info: part of spreading distribution
           name: name of it
        Returns:
           dist: dist hash
        """
        dist = {}
        info = info.rstrip("_")
        splitted = info.split("_")
        if name == Trace.SPROB:
           assert len(splitted) <= 2
           start = float(splitted[0])
           if len(splitted) == 2:
              end = float(splitted[1])
              dist[name] = {"start":start, "end":end}
           else:
              dist[name] = start  
        else:
           assert splitted[0] in ["expo","rayleigh","powerlaw","weibull","rayleigh"]
           if splitted[0] in ["expo","rayleigh","powerlaw"]:
              assert len(splitted) in [2,3]
              if len(splitted) == 2:
                 dist[name] = (splitted[0],(float(splitted[1]),))
              else:
                 dist[name] = {"dist": splitted[0], "start": (float(splitted[1]),), "end":(float(splitted[2]),) }   
           elif splitted[0] in ["weibull","lognormal"]:
              assert len(splitted) in [3,5]
              if len(splitted) == 3:
                 start,end = [float(item) for item in splitted[1:]]  
                 dist[name] = (splitted[0],(start,end))
              else:
                 start1,start2,end1,end2 = [float(item) for item in splitted[1:]] 
                 dist[name] = {"dist": splitted[0], "start": (start1,start2), "end":(end1,end2) }   
        return dist
    
    @classmethod
    def parseTransition(self,info,name):
        """ parse given transition string to distribution
        Args:
           info: part of spreading distribution
           name: name of it
        Returns:
           dist: dist hash
        """
        dist = {}
        info = info.rstrip("_")
        splitted = info.split("_")
        if name == Trace.SPROB:
           assert len(splitted) <= 2
           start = float(splitted[0])
           if len(splitted) == 2:
              end = float(splitted[1])
              dist[name] = {"start":start, "end":end}
           else:
              dist[name] = start  
        else:
           assert splitted[0] in ["expo","rayleigh","powerlaw","weibull","rayleigh"]
           if splitted[0] in ["expo","rayleigh","powerlaw"]:
              assert len(splitted) in [2,3]
              if len(splitted) == 2:
                 dist[name] = (float(splitted[1]),)
              else:
                 dist[name] = {"dist": splitted[0], "start": (float(splitted[1]),), "end":(float(splitted[2]),) }   
           elif splitted[0] in ["weibull","lognormal"]:
              assert len(splitted) in [3,5]
              if len(splitted) == 3:
                 start,end = [float(item) for item in splitted[1:]]  
                 dist[name] = (start,end)
              else:
                 start1,start2,end1,end2 = [float(item) for item in splitted[1:]] 
                 dist[name] = {"dist": splitted[0], "start": (start1,start2), "end":(end1,end2) }    
        return dist
    
    @classmethod
    def folder2SpreadInfo(self,folder):
        """converts given foldername to spread info trace folder may also include start and end
        Args:
          folder:
        Returns:
          smodel:
          prob:
          distinfo: 
        """
        distinfo = {}
        prob,smodel,info = folder.split("-")
        items =  [spread + "_" for spread in [Trace.S2I,Trace.S2E,Trace.I2R,Trace.E2I,Trace.I2S,Trace.SPROB] if info.find(spread) != -1]
        for item1 in items:
            infodict = {}
            for item2 in items:
                if item1 == item2:
                   continue
                mystr = info.split(item1)[1].split(item2)[0]
                infodict[mystr] = len(mystr)
            dist = self.parseTransitionWithDist(min(infodict.items(), key=lambda x: x[1])[0],item1[0:-1])
            key = dist.keys()[0]
            distinfo[key] = dist[key]
        return smodel,prob,distinfo
    
    @classmethod
    def getSpreadFolder(self,smodel,dists,dist):
       """returns foldername 
       Args:
         smodel: spreading model
         dists: distributions
         dist: continuous/discrete
       Returns:
         sfolder: spreading foldername
       """
       arr = []
       for key in dists.keys():
           if key in [self.I2R, self.E2I, self.S2I, self.I2S, self.S2E]:
              arr.append(key)
              if type(dists[key]) == dict:
                 arr.append(str(dists[key]["dist"]))
                 arr.extend([str(item2) for item in ["start","end"] for item2 in list(dists[key][item])])
              else:
                 arr.extend([str(dists[key][0])] + [str(item) for item in list(dists[key][1])])
       if dists.has_key(self.SPROB):
          arr.append(self.SPROB)
          if type(dists[key]) == dict:
             arr.extend([str(dists[self.SPROB]["start"]),str(dists[self.SPROB]["end"])])
          else:
             arr.append(str(dists[self.SPROB]))
       return "{0}-{1}-{2}".format(dist,smodel,"_".join(arr))

    @classmethod
    def getRandomTraceFile(self,tracefolder,fraccons):
       """ returns random trace file under tracefolder having at least minicount infected nodes
       Args:
          tracefolder:
          fraccons:
       Returns:
          tracefile:
       """
       tracefiles = myutil.listfiles(tracefolder)
       random.shuffle(tracefiles)
       for tracefile in tracefiles:
           tracepath = "{0}/{1}".format(tracefolder,tracefile) 
           from InputOutput import InputOutput #!!!!!
           trace = InputOutput.readPlainTrace(tracepath)
           inodes = [node for node in trace.keys() if trace[node].has_key(Trace.INFECTED)]
           if len(inodes) >= fraccons[Trace.INFECTED]["min"] and len(inodes) <= fraccons[Trace.INFECTED]["max"]:
              return tracefile
       print "None of traces under {0} has between {1} and {2} infected nodes".format(tracefolder,fraccons[Trace.INFECTED]["min"],fraccons[Trace.INFECTED]["max"])
       print "EXITING !!!"
       exit(1)
       
    @classmethod
    def getMaxTraceTime(self,trace):
        """gets maximum trace time
        Args:
           trace: trace data
        Returns:
           maxtime: maxtime
        """
        return max([trace[node][state] for node in trace.keys() for state in trace[node].keys()])

    @classmethod 
    def trace2Snapshot(self,trace,obstime,smodel,G):
       """returns snapshot at timefrac given trace data
       Args:
         trace: trace data
         obstime: given time
         smodel: spreading model
         G: graph
       Returns:
         curstate: returns current state out of given trace
       """
       maxtime = Trace.getMaxTraceTime(trace) 
       curstate = {self.INFECTED: set(), self.SUSCEPTIBLE: set([node for node in G.nodes() if not trace.has_key(node)]), self.RECOVERED: set(), self.EXPOSED: set()}
       for node in G.nodes():
           if not trace.has_key(node):
              curstate[self.SUSCEPTIBLE].add(node)
              continue
           etime, itime, rtime = maxtime*100000, maxtime*100000, maxtime*100000
           if trace[node].has_key(self.EXPOSED):
              etime = trace[node][self.EXPOSED]
           if trace[node].has_key(self.INFECTED):
              itime = trace[node][self.INFECTED]
           if trace[node].has_key(self.RECOVERED):
              rtime = trace[node][self.RECOVERED]
           if smodel in ["si", "sis","sir"] and obstime < itime:
              curstate[self.SUSCEPTIBLE].add(node)
           elif smodel == "seir" and obstime < etime and obstime < itime:
              curstate[self.SUSCEPTIBLE].add(node)
           elif obstime >= etime and obstime < itime:
              curstate[self.EXPOSED].add(node)
           elif obstime >= itime and obstime < rtime:
              curstate[self.INFECTED].add(node)
           elif obstime >= rtime:
              curstate[self.RECOVERED].add(node)
       return curstate

    @classmethod 
    def addStateChangeNoise(self,curstate,noisedeg,smodel):
       """add state change noise to a given snapshot
         Args:
           curstate: current state
           noisedeg: noise degree
           smodel: spreading model
         Returns:
           curstate: returns new noise added current state
       """
       statemap = {0: self.SUSCEPTIBLE, 1: self.INFECTED}
       if smodel in ["sir","seir"]:
          statemap[2] = self.RECOVERED
       if smodel == "seir":
          statemap[3] = self.EXPOSED
       changenodes = {state : set(node for node in curstate[state] if random.random() <= noisedeg) for state in curstate.keys()}
       for state in changenodes.keys():
           for node in changenodes[state]:
               rnum = random.randint(0,len(statemap.keys())-1)
               newstate,oldstate = statemap[rnum], state
               curstate[oldstate].remove(node)
               curstate[newstate].add(node)
       return curstate

    @classmethod       
    def addNoise(self,curstate,noisedeg,noisetype,smodel):
       """noises a given snapshot
         Args:
           curstate: current state
           noisedeg: noise degree
           noisetype: type of noise
           smodel: spreading model
         Returns:
           curstate: returns new noise added current state
       """
       assert noisetype in ["Partial", "StateChange"]
       if noisedeg == 0.0:
          return curstate
       if noisetype == "Partial":
          print "partial noise not impl yet!!"
       elif noisetype == "StateChange":
          noisemethod = "add{0}Noise".format(noisetype)
          method = getattr(self,noisemethod)
          return method(curstate,noisedeg,smodel)

    @classmethod
    def addStateChangeNoise2Trace(self,trace,noise,smodel,G):
        """adds State change noise to trace
        Args:
           trace:
           noise:
           smodel:
           G:
        Returns:
           newtrace:
        """
        maxtime = max([trace[node][state] for node in trace.keys() for state in trace[node].keys()])
        statemap = {"si":[Trace.INFECTED],"sir":[Trace.INFECTED,Trace.RECOVERED],"seir":[Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]}
        INFTIME,newtrace,startnodes = {}, {}, []
        if smodel == "si":
           statepairs, nonstatepairs = [], []
           INFTIME[Trace.INFECTED] =  maxtime + 1
        elif smodel == "sir":
           statepairs = [(Trace.INFECTED, Trace.RECOVERED)]
           nonstatepairs = [(Trace.RECOVERED, Trace.INFECTED)]
           INFTIME ={Trace.INFECTED:maxtime+1,Trace.RECOVERED:maxtime + 2}
        elif smodel == "seir":
           statepairs = [(Trace.EXPOSED, Trace.INFECTED),(Trace.INFECTED, Trace.RECOVERED),(Trace.EXPOSED, Trace.RECOVERED)]
           nonstatepairs = [(Trace.RECOVERED, Trace.INFECTED)]
           INFTIME ={Trace.EXPOSED:maxtime+1, Trace.INFECTED:maxtime+2,Trace.RECOVERED:maxtime + 3}
           startnodes = [node for node in trace.keys() if trace[node].has_key(Trace.INFECTED) and trace[node][Trace.INFECTED] == 0]
        for node in G.nodes():
            if node in startnodes:
               newtrace[node] = deepcopy(trace[node])  
               continue
            newtrace[node] = {}
            while True:
                for state in statemap[smodel]:
                    if not trace.has_key(node) or not trace[node].has_key(state):
                       curtime = INFTIME[state]
                    else:
                       curtime = trace[node][state]
                    randfrac = random.random() 
                    if randfrac <= noise:
                       extra = int(round(randfrac * maxtime))
                       if random.random() <= 0.5:
                          extra *= -1
                       newtime = max(0, extra+curtime)
                    else:
                       newtime = curtime
                    if newtime <= maxtime:   
                       newtrace.setdefault(node,{})
                       newtrace[node][state] = newtime
                flag = True  
                for state1,state2 in statepairs:
                    if newtrace[node].has_key(state1) and newtrace[node].has_key(state2):
                       if newtrace[node][state1] >= newtrace[node][state2]:
                          flag = False
                          break
                for state1,state2 in nonstatepairs:       
                    if newtrace[node].has_key(state1) and not newtrace[node].has_key(state2):
                       flag = False
                       break
                if flag:
                   break
            if len(newtrace[node].keys()) == 0:
               del newtrace[node] 
        return newtrace
    
    @classmethod  
    def modifyTraces(self,traces,samplerate,noiseinfo,smodel):
        """ modifies traces (subsampling noise addition etc) (always add probabilistic noise)
        Args:
           traces: traces 
           samplerate: sampling rate
           noiseinfo: noise information
           smodel: spreading model
        Returns:
           newtraces: new traces
        """
        noiseshape,noise = noiseinfo
        newtraces = []
        for trace in traces:
            sampledtrace = self.sampleTrace(trace,smodel,samplerate)
            if noise != 0.0:
               newtraces.append(self.makeTracePartial(sampledtrace,noiseinfo,smodel)) 
            else:
               newtraces.append(sampledtrace)
        return newtraces
 
    @classmethod
    def sampleTrace(self,trace,smodel,samplerate):
       """subsamples trace data
       Args:
         trace: trace data
         smodel: spreading model
         samplerate: samplerate(0 is perfect case)
       Returns:
         sampledtrace:
       """
       assert smodel in ["sir","si","seir"]
       if samplerate == 0:
          return trace

       startnode,starttime = [(node,min(trace[node].values())) for node in trace.keys() if len(set(trace[node].values()).intersection(set([0,0.0]))) != 0 ][0]
       subtrace = {node: {state: int(math.ceil(float(trace[node][state])/samplerate)) * samplerate for state in trace[node].keys()} for node in trace.keys() }
       if smodel == "sir":               
          for node in subtrace.keys():
              if not subtrace[node].has_key(self.RECOVERED):
                 continue
              assert subtrace[node][self.INFECTED] <= subtrace[node][self.RECOVERED]
              if subtrace[node][self.INFECTED] == subtrace[node][self.RECOVERED]:
                 del subtrace[node][self.INFECTED]
       elif smodel == "seir":               
          for node in subtrace.keys():
              if not subtrace[node].has_key(self.INFECTED) or not subtrace[node].has_key(self.RECOVERED):
                 continue
              assert subtrace[node][self.EXPOSED] <= subtrace[node][self.INFECTED] and subtrace[node][self.INFECTED] <= subtrace[node][self.RECOVERED]
              if subtrace[node][self.INFECTED] == subtrace[node][self.EXPOSED]:
                 del subtrace[node][self.EXPOSED]  
              if subtrace[node][self.INFECTED] == subtrace[node][self.RECOVERED]:
                 del subtrace[node][self.INFECTED]
       subtrace[startnode][Trace.INFECTED] = 0
       return subtrace

    @classmethod
    def assignNodeStates(self,trace,smodel):
       """converts tracedata into format for noisy case
       Args:
          trace: trace data
          smodel: spreading model
       Returns:
          returns new data structure for trace
       """
       assert smodel in ["si", "sir", "seir"]
       nodestates = {}
       alltimes = sorted(list(set(trace[node][state] for node in trace.keys() for state in trace[node].keys())))
       if smodel == "si":
          nodestates = {node: {time: {self.SUSCEPTIBLE:1, self.INFECTED:0} for time in alltimes} for node in trace.keys()}
          for node in trace.keys():
              for time in alltimes:
                  if time >= trace[node][self.INFECTED]:
                     nodestates[node][time] = {self.SUSCEPTIBLE:0, self.INFECTED:1}
       elif smodel == "sir":       
          for node in trace.keys():
              nodestates[node] = {}
              itime,rtime = 10000000000, 10000000004 #very large numbers
              if trace[node].has_key(self.INFECTED):
                 itime = trace[node][self.INFECTED]
              if trace[node].has_key(self.RECOVERED):
                 rtime = trace[node][self.RECOVERED]
              for time in alltimes:
                  nodestates[node][time] = {self.SUSCEPTIBLE:1, self.INFECTED:0, self.RECOVERED:0}
                  if time >= itime and time < rtime:
                     nodestates[node][time] = {self.SUSCEPTIBLE:0, self.INFECTED:1, self.RECOVERED:0} 
                  elif time >= rtime:
                     nodestates[node][time] = {self.SUSCEPTIBLE:0, self.INFECTED:0, self.RECOVERED:1} 
       elif smodel == "seir":
          for node in trace.keys():
              nodestates[node] = {}
              etime, itime, rtime = 10000000000, 10000000002, 10000000004 #very large number
              if trace[node].has_key(self.EXPOSED):
                 etime = trace[node][self.EXPOSED]
              if trace[node].has_key(self.INFECTED):
                 itime = trace[node][self.INFECTED]
              if trace[node].has_key(self.RECOVERED):
                 rtime = trace[node][self.RECOVERED]
              for time in alltimes:
                  nodestates[node][time] = {self.SUSCEPTIBLE:1, self.EXPOSED:0, self.INFECTED:0, self.RECOVERED:0} 
                  if time >= etime and time < itime:
                     nodestates[node][time] = {self.SUSCEPTIBLE:0, self.EXPOSED:1, self.INFECTED:0, self.RECOVERED:0}
                  elif time >= itime and time < rtime:
                     nodestates[node][time] = {self.SUSCEPTIBLE:0, self.EXPOSED:0, self.INFECTED:1, self.RECOVERED:0}
                  elif time >= rtime:
                     nodestates[node][time] = {self.SUSCEPTIBLE:0, self.EXPOSED:0, self.INFECTED:0, self.RECOVERED:1}
       return nodestates

    @classmethod
    def makeTracePartial(self,trace,noiseinfo,smodel):
       """makes given data noisy
       Args:
          trace: trace data
          noiseinfo: noise 
          smodel: spreading model
       Returns:
          returns noisy trace data
       """
       noiseshape,noise = noiseinfo
       if noiseshape == "add":
          return self.makeTracePartialAdditive(trace,noise,smodel)
       else:
          return self.makeTracePartialMultiplicative(trace,noiseinfo,smodel)

    @classmethod
    def genKernelMatrix(self,type,maxtime,noise=0.1):
        """generates kernel matrix
        Args:
           type: type of kernel
           maxtime: kernel leng will be 2*maxtime
           rho2 = rho2 for dispersion
        Returns:
           kernmat: kernel matrix
        """
        assert type in ["gaussian"]
        if type == "gaussian":
           rho2 = 100 * noise**2
           kernmat = np.zeros((2*maxtime+1,2*maxtime+1),dtype = np.float64)
           for time1 in xrange(np.shape(kernmat)[0]):
               for time2 in xrange(np.shape(kernmat)[1]):
                   kernmat[time1,time2] = math.exp((-0.5 * float(time1 - time2)**2) /rho2)
        return kernmat
       
    @classmethod
    def makeTracePartialMultiplicative(self,trace,noiseinfo,smodel):
       """makes given data noisy by matrix multiplicative noise
       Args:
          trace: trace data
          noiseinfo: noise info about noise and shape 
          smodel: spreading model
       Returns:
          returns noisy trace data
       """
       assert smodel in ["si"]
       noiseshape, noise = noiseinfo
       assert noiseshape in ["uniform","gauss","multigauss"]
       nodestates = self.assignNodeStates(trace,smodel)
       sortedtimes = sorted(list(set([time for node in nodestates.keys() for time in nodestates[node].keys()])))
       if noiseshape == "uniform":
          shapemat = np.zeros((len(sortedtimes),len(sortedtimes)),dtype=np.float64)
          univec = [random.random()]
          for index in xrange(len(sortedtimes)-1):
              if random.random() <= 1.0 - noise:
                 univec.append(univec[-1])
              else:    
                 univec.append(random.random())
          univec = np.array(univec)
          for index1 in xrange(len(sortedtimes)):
              shapemat[index1,:] = univec
          if smodel == "si":
             for node in nodestates.keys():
                 pervec = np.zeros((len(sortedtimes),),dtype=np.float64)
                 if trace[node].has_key(self.INFECTED):
                    tindex = sortedtimes.index(trace[node][self.INFECTED])  
                    pervec[tindex] = 1.0
                 endvec = np.dot(shapemat,pervec)
                 for tindex in xrange(len(sortedtimes)):
                     time = sortedtimes[tindex]
                     nodestates[node][time] = {self.SUSCEPTIBLE: 1.0 - endvec[tindex], self.INFECTED: endvec[tindex] }
             for node in nodestates.keys():
                 for time in nodestates[node].keys():
                     assert nodestates[node][time][self.INFECTED] <= 1.0
       elif noiseshape == "gauss":
          kernmat = self.genKernelMatrix("gaussian",sortedtimes[-1],noise)
          shapemat = np.zeros((len(sortedtimes),len(sortedtimes)),dtype=np.float64)
          if smodel == "si":
             for node in nodestates.keys():
                 ivec = np.zeros((len(sortedtimes),),dtype=np.float64)
                 flag = False
                 if trace[node].has_key(self.INFECTED):
                    tindex = sortedtimes.index(trace[node][self.INFECTED])  
                    ivec[tindex] = 1.0
                    flag = True
                 if flag:
                    startpos = (-1.0 * trace[node][self.INFECTED]) + (np.shape(kernmat)[0] / 2)
                    for index1 in xrange(len(sortedtimes)):
                        pos1 = sortedtimes[index1] - startpos
                        for index2 in xrange(len(sortedtimes)):
                            pos2 = sortedtimes[index2] - startpos
                            shapemat[index1,index2] = kernmat[pos1,pos2]
                    maxval = shapemat.max()
                    shapemat /= float(maxval)
                    ivec = np.multiply(shapemat,ivec)
                 for tindex in xrange(len(sortedtimes)):
                    time = sortedtimes[tindex]
                    nodestates[node][time] = {self.SUSCEPTIBLE: 1.0 - ivec[tindex], self.INFECTED: ivec[tindex]}
       elif noiseshape == "multigauss":
          kernmat = self.genKernelMatrix("gaussian",sortedtimes[-1])
          shapemat = np.zeros((len(sortedtimes),len(sortedtimes)),dtype=np.float64)
          print "not implemeted yet!!"
          exit(1)
       return nodestates
        
    @classmethod
    def makeTracePartialAdditive(self,trace,noise,smodel):
       """makes given data noisy additive noise
       Args:
          trace: trace data
          noise: noise 
          smodel: spreading model
       Returns:
          returns noisy trace data
       """
       assert smodel in ["si", "sir", "seir"]
       nodestates = self.assignNodeStates(trace,smodel)
       alltimes = sorted(list(set([time for node in nodestates.keys() for time in nodestates[node].keys()])))
       if smodel == "si":
          for time in alltimes:
              for node in nodestates.keys():
                  probs = nodestates[node][time]
                  if random.random() < (1.0-noise)**2:
                     continue
                  corrupt = random.uniform(0,noise)
                  if probs[self.SUSCEPTIBLE] == 1:
                     nodestates[node][time] = {self.SUSCEPTIBLE: 1.0 - corrupt, self.INFECTED: corrupt }
                  elif probs[self.INFECTED] == 1:
                     nodestates[node][time] = {self.SUSCEPTIBLE: corrupt, self.INFECTED: 1.0 - corrupt }
       elif smodel == "sir":
          for time in alltimes:
              for node in nodestates.keys():
                  probs = nodestates[node][time]
                  if random.random() < (1.0-noise)**2:
                     continue
                  corrupt = random.uniform(0,noise)   
                  if probs[self.SUSCEPTIBLE] == 1:
                     nodestates[node][time] = {self.SUSCEPTIBLE: 1.0 - corrupt, self.INFECTED: corrupt, self.RECOVERED: 0.0}
                  elif probs[self.INFECTED] == 1:
                     if random.random() <= 0.5:
                        nodestates[node][time] = {self.SUSCEPTIBLE: corrupt, self.INFECTED: 1.0 - corrupt, self.RECOVERED: 0.0 }
                     else:
                        nodestates[node][time]= {self.RECOVERED: corrupt, self.INFECTED: 1.0 - corrupt, self.SUSCEPTIBLE: 0.0 }
                  elif probs[self.RECOVERED] == 1:
                     nodestates[node][time] = {self.RECOVERED: 1.0 - corrupt, self.INFECTED: corrupt, self.SUSCEPTIBLE: 0.0 }
       elif smodel == "seir":
          for time in alltimes:
              for node in nodestates.keys():
                  probs = nodestates[node][time]
                  if random.random() < (1.0-noise)**2:
                     continue
                  corrupt = random.uniform(0,noise)   
                  if probs[self.SUSCEPTIBLE] == 1:
                     nodestates[node][time] = {self.SUSCEPTIBLE: 1.0 - corrupt, self.EXPOSED: corrupt, self.INFECTED: 0.0, self.RECOVERED: 0.0 }
                  elif probs[self.EXPOSED] == 1:
                     if random.random() <= 0.5:
                        nodestates[node][time] = {self.SUSCEPTIBLE: corrupt, self.EXPOSED: 1.0 - corrupt, self.INFECTED: 0.0, self.RECOVERED: 0.0 }
                     else:
                        nodestates[node][time]= {self.INFECTED: corrupt, self.EXPOSED: 1.0 - corrupt, self.SUSCEPTIBLE: 0.0, self.RECOVERED: 0.0 }
                  elif probs[self.INFECTED] == 1:
                     if random.random() <= 0.5:
                        nodestates[node][time] = {self.EXPOSED: corrupt, self.INFECTED: 1.0 - corrupt, self.SUSCEPTIBLE: 0.0, self.RECOVERED: 0.0 }
                     else:
                        nodestates[node][time] = {self.RECOVERED: corrupt, self.INFECTED: 1.0 - corrupt, self.SUSCEPTIBLE: 0.0, self.EXPOSED: 0.0 }
                  elif probs[self.RECOVERED] == 1:
                     nodestates[node][time] = {self.RECOVERED: 1.0 - corrupt, self.INFECTED: corrupt, self.SUSCEPTIBLE: 0.0, self.EXPOSED: 0.0 }
       return nodestates

