##History Ensemble
import os
import sys
sys.path.append("./lib")
import myutilities as myutil
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from copy import deepcopy


class DiscreteNoneEnsemble():
     """ Ensemble Solution Class for Discrete No Bound(used also for bound case)
     """
     def __init__(self):
        return

     @classmethod
     def getAbsDif(self,eststate,curstate):
        """difference between eststate and curstate
        Args:
          eststate:
          curstate:
        Returns:
          dif: sum difference
        """
        dif = 0.0
        for item in [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]:
            if eststate.has_key(item) and curstate.has_key(item):
               dif += len(eststate[item].difference(curstate[item]))
               dif += len(curstate[item].difference(eststate[item]))
            elif eststate.has_key(item) and not curstate.has_key(item):
               dif += len(eststate[item])
            elif not eststate.has_key(item) and curstate.has_key(item):
               dif += len(curstate[item])
        return dif/2

     @classmethod
     def getSubObj(self,startnodes,curstates,smodel,obstimes,G):
        """returns objective value submodular
        Args:
          startnodes: start nodes
          curstates: current state
          smodel: spreading model
          obstimes: 
          G: 
        Returns:
          avgdif: average dif
        """
        assert len(obstimes) == len(curstates)
        weights = [1.0]
        for index in xrange(len(obstimes)-1):
            weights.append(weights[index] * 0.9)
        ITERCOUNT,dif = 5, 0.0
        for index in xrange(ITERCOUNT):
            tracemethod = "gen{0}Trace".format(smodel.capitalize())
            method = getattr(DisTrace,tracemethod)
            (stimes,etimes,itimes,rtimes) = method(G,startnodes,"static")
            trace = Trace.convertTimes2Trace(smodel,stimes,etimes,itimes,rtimes)
            for index in xrange(len(obstimes)):
                obstime = obstimes[index]
                curstate = curstates[index]
                eststate = Trace.trace2Snapshot(trace,obstime,smodel,G)
                dif += self.getAbsDif(eststate,curstate)
        return dif / (float(ITERCOUNT) * sum(weights))

     @classmethod    
     def GreedyAddition(self,nodes,obstimes,curstates,smodel,G):
        """ greedily estimates the subset which will best describe the the observed trace
        Args:
           nodes: nodes at time: time
           obstimes: obs times
           curstates: current states
           smodel: spreading model
           G:
        Returns:
           sol: solution 
           objval: objective value
        """
        objval, curnodes, sol = 100000000.0, set(nodes), set()
        while len(curnodes) >= 0:
           decr = {}  
           for node in curnodes:
               tempsol = sol.union(set([node]))
               decrement = objval - self.getSubObj(tempsol,curstates,smodel,obstimes,G)
               print "{0}: dec {1}".format(tempsol,decrement)
               if decrement > 0:
                  decr[node] = decrement
           if len(decr.keys()) == 0:
              break
           node = max(decr.items(), key=lambda x: x[1])[0]
           sol.add(node)
           curnodes.remove(node)
           objval -= decr[node]
        return sol,objval

     @classmethod
     def getSpreaderEnsembleSolution(self,allsol,curstates,obstimes,G,smodel,inter):
        """returns ensemble solution
          Args:
             allsol: solutions
             curstates: current state knowledge
             obstimes: obstimes
             G: Graph
             smodel: spreading model
             inter:
          Returns:
             bestnodes: solution
             besttime
        """
        newcurstates = deepcopy(curstates)
        time2nodes = {}
        for hist in allsol:
            #mintime = None
            #for time in sorted(hist.keys()):
            #    flag = True
            #    for state in curstates[0].keys():
            #        if len(set(curstates[0][state]).symmetric_difference(hist[time][state])) != 0:
            #           flag = False
            #           break
            #    if flag:
            #       mintime = time
            #       break
            #assert mintime != None 
            mintime = min(hist.keys())
            time2nodes.setdefault(mintime,set())
            time2nodes[mintime] |= set(hist[mintime][Trace.INFECTED])
        if inter == "bound":  
           assert time2nodes.keys() == [0] and type(curstates) == dict
           newcurstates,obstimes = [],[]
           for time in sorted(curstates.keys()):
               obstimes.append(time)
               newcurstates.append(deepcopy(curstates[time]))
           assert type(newcurstates) == list        
        bestscore, besttime, bestnodes = 10000000000.0, -1, set()
        for curtime in time2nodes.keys():
            print "curtime:",curtime
            if inter == None:
               sentobstimes = [obstime-mintime for obstime in obstimes] 
            curnodes = set(time2nodes[curtime])  
            retnodes,score = self.GreedyAddition(curnodes,sentobstimes,newcurstates,smodel,G)
            if score <= bestscore:
               bestscore,besttime,bestnodes = score,curtime,retnodes
            continue  #if this is not here, loop will stop after first iteration!!!  
	    assert len(bestnodes) > 0 
        return {Trace.INFECTED: set(bestnodes), Trace.SUSCEPTIBLE:set(), Trace.EXPOSED:set(), Trace.RECOVERED:set()},besttime
     
     @classmethod
     def getEnsembleSolution(self,allsol,infermode,curstates,obstimes,G,smodel,inter=None):
         """returns ensemble solution
           Args:
             allsol: solutions
             infermode: infermode
             curstates: current state knowledge
             obstimes:
             G: Graph
             smodel: spreading model
             inter:
           Returns:
             bestnodes: solution
             besttime
         """
         assert inter in [None,"bound"]
         assert infermode == "Spreader"
         if infermode == "Spreader":
            return self.getSpreaderEnsembleSolution(allsol,curstates,obstimes,G,smodel,inter)
         
