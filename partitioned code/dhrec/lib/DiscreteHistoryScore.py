import networkx as nx
import numpy as np
import scipy as sp
import scipy.stats
import random
from copy import deepcopy
import math
import itertools
import sys
import os
import myutilities as myutil
from DiscreteTrace import DiscreteTrace as Trace
from InputOutput import InputOutput
import checkUtil

class DiscreteHistoryScore():
    """history inference score class for discrete distribution
    """
    def __init__(self):
        return

    @staticmethod
    def trace2History(trace,evol,smodel,giventimes):
        """converts trace to dict format for score estimation(negates times and deletes giventimes)
        Args:
          trace: trace
          evol: static/dynamic
          smodel: spreading model
          giventimes: given trace times
        Returns:
          sol: dict format for score estimation 
        """
        trstates = set(trace[node][state] for node in trace.keys() for state in trace[node].keys())
        sol,modsol = {},{}
        for node in trace.keys():
            for state in trace[node].keys():
                sol.setdefault(trace[node][state],{})
                sol[trace[node][state]].setdefault(state,set()).add(node)
        maxtime = max(giventimes)
        for time in sol.keys():
            if time > maxtime:
               del sol[time]
        for time in xrange(0,maxtime+1):
            if not sol.has_key(time):
               sol[time] = {state:set() for state in trstates}
        curtimes = sol.keys()       
        for time in curtimes:
            modsol[time - maxtime] = deepcopy(sol[time])
        return modsol

    @staticmethod
    def hist2Times(hist,smodel):
        """converts hist data to state times
        Args:
           hist:
           smodel:
        Returns:
           stimes,etimes,itimes,rtimes:
        """
        times = sorted([time for time in hist.keys() if time <= obstime])
        allnodes = set(node for time in hist.keys() for state in hist[time].keys() for node in hist[time][state]) 
        statetimes = {Trace.SUSCEPTIBLE:{},Trace.EXPOSED:{},Trace.INFECTED:{},Trace.RECOVERED:{}}
        for node in allnodes:
            for state in [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]:
                foundtime = None 
                for time in times:
                    if hist[time].has_key(state) and node in hist[time][state]:
                       foundtime = time
                       break
                if foundtime != None:
                   statetimes[state][node] = foundtime     
        return statetimes[Trace.SUSCEPTIBLE],statetimes[Trace.EXPOSED],statetimes[Trace.INFECTED],statetimes[Trace.RECOVERED]

    @staticmethod
    def hist2SnapshotWithGraph(hist,smodel,G,obstime):
        """converts hist data to snapshot at obstime with graph given
        Args:
           hist:
           smodel:
           G:
           obstime:
        Returns:
           curstate:
        """
        curstate = DiscreteHistoryScore.hist2Snapshot(hist,smodel,obstime)
        seennodes = set(node for state in curstate.keys() for node in curstate[state])
        curstate[Trace.SUSCEPTIBLE] |= set(G.nodes()).difference(seennodes)
        assert sum([len(curstate[state]) for state in curstate.keys()]) == G.number_of_nodes() 
        return curstate
            
    @staticmethod
    def hist2Snapshot(hist,smodel,obstime):
        """converts hist data to snapshot at obstime
        Args:
           hist:
           smodel:
           obstime:
        Returns:
           curstate:
        """
        times = sorted([time for time in hist.keys() if time <= obstime])
        allnodes = set(node for time in hist.keys() for state in hist[time].keys() for node in hist[time][state]) 
        curstate = {Trace.SUSCEPTIBLE:set(),Trace.EXPOSED:set(),Trace.INFECTED:set(),Trace.RECOVERED:set()}
        for node in allnodes:
            foundstate = Trace.SUSCEPTIBLE
            for state in [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]:
                for time in times:
                    if hist[time].has_key(state) and node in hist[time][state]:
                       foundstate = state
                       break
            curstate[foundstate].add(node)        
        return curstate

    @staticmethod
    def getGraphBasedRandScore(dict1,dict2,items):
        """gets graph based rand score
        Args:
           dict1:
           dict2:
           items:      
        Returns:
           randscore:
        """
        statemap1 = {node:state for state in dict1.keys() for node in dict1[state]}
        statemap2 = {node:state for state in dict2.keys() for node in dict2[state]}
        a,b,c,d = 0,0,0,0
        for index1 in xrange(len(items)):
            node1 = items[index1]
            for index2 in xrange(index1+1,len(items)):
                node2 = items[index2]     
                if statemap1[node1] == statemap1[node2] and statemap2[node1] == statemap2[node2]: 
                   a += 1 
                elif statemap1[node1] != statemap1[node2] and statemap2[node1] != statemap2[node2]: 
                   b += 1 
                elif statemap1[node1] != statemap1[node2] and statemap2[node1] == statemap2[node2]: 
                   c += 1 
                elif statemap1[node1] == statemap1[node2] and statemap2[node1] != statemap2[node2]: 
                   d += 1 
        assert a+b+c+d == (len(items) * (len(items)-1))/2             
        return float(a+b) / (a+b+c+d)                                  

    @staticmethod 
    def getGraphBasedJaccardSim(dict1,dict2,states):
        """gets graph based jaccard sim
        Args:
           dict1:
           dict2: 
           states:   
        Returns:
           avgjaccardscore:
        """
        curdif = 0.0
        for state in states:
            curset1,curset2 = set(), set()
            if dict1.has_key(state):
               curset1 = set(dict1[state])  
            if dict2.has_key(state):
               curset2 = set(dict2[state])   
            interset = curset1.intersection(curset2)
            unionset = curset1.union(curset2)
            if len(unionset) == 0:
               curdif += 1.0
            else: 
               curdif += float(len(interset))/len(unionset)
        return curdif / float(len(states)) 
    
    @staticmethod 
    def getPartGraphScore(truth,sol,G,iternum,scoretype,smodel,evol,obstimes,inter):
        """history prediction score also using graph structure #???
        Args:
          truth: truth dict
          sol: our solution dict
          G: Graph
          iternum: iteration count
          scoretype: scoretype
          smodel: spreading model
          evol: static/dynamic
          obstimes:
          inter:
        Returns:
          score: avg graphwise closeness score 
        """
        for obstime in obstimes:
            assert obstime >= 0
        assert scoretype in ["jaccardsim","rand"]    
        if smodel == "si":
           states = [Trace.SUSCEPTIBLE,Trace.INFECTED]
        elif smodel == "si":
           states = [Trace.SUSCEPTIBLE,Trace.INFECTED,Trace.RECOVERED]
        elif smodel == "seir":
           states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]  
        truthsnap = {} #truth in terms of snapshots
        for time in sorted(truth.keys()):
            truthsnap[time] = DiscreteHistoryScore.hist2SnapshotWithGraph(truth,smodel,G,time)
            assert sum([len(truthsnap[time][state]) for state in truthsnap[time].keys()]) == G.number_of_nodes()
        if inter in ["None",None]:
           add = -1 * min(truthsnap.keys())
           assert min(truthsnap.keys()) < 0 and max(obstimes) + min(truthsnap.keys()) == 0
        elif inter == "bound":
           add = 0 
           assert min(truthsnap.keys()) >= 0
        knowns = [deepcopy(truthsnap[obstime - add]) for obstime in obstimes]
        for time in sorted(sol.keys()):
            if len(sol[time][Trace.INFECTED]) != 0:
               initnodes = set(sol[time][Trace.INFECTED])   
               break
        tracemethod = "gen{0}Trace".format(smodel.capitalize())
        method = getattr(Trace,tracemethod)
        totdif = 0.0
        for index in xrange(iternum):
            (stimes,etimes,itimes,rtimes) = method(G,initnodes,evol)
            trace = Trace.convertTimes2Trace(smodel,stimes,etimes,itimes,rtimes)
            cursolstate = []
            maxtime = max([trace[node][state] for node in trace.keys() for state in trace[node].keys()])
            for obstime in obstimes:
                if obstime >= maxtime:
                   solstate = Trace.trace2Snapshot(trace,maxtime,smodel,G) 
                else:
                   solstate = Trace.trace2Snapshot(trace,obstime,smodel,G)
                cursolstate.append(deepcopy(solstate))
            for index in xrange(len(cursolstate)):
                if scoretype == "rand":
                   totdif += DiscreteHistoryScore.getGraphBasedRandScore(cursolstate[index],knowns[index],G.nodes()) 
                elif scoretype == "jaccardsim":
                   totdif += DiscreteHistoryScore.getGraphBasedJaccardSim(cursolstate[index],knowns[index],states)  
        return totdif / (float(iternum) * len(cursolstate))              

    @staticmethod  
    def getMatchingScore(truth,sol,G,usestates):
        """estimates average matching score
        Args:
           truth:
           sol:
           G:
           usestates:
        Returns:
           mscore:
        """
        assert len(usestates) == 1
        state = usestates[0]
        diam = -10.0
        for node1 in G.nodes():
            for node2 in G.nodes():
                try:
                   dist = nx.shortest_path_length(G,source=node1,target=node2)
                except:
                   pass
                if dist > diam:
                   diam = dist
        mintruth = min([time for time in truth.keys() if truth[time].has_key(state)])
        minsol = min([time for time in sol.keys() if sol[time].has_key(state)])
        matchG = nx.Graph()
        for node1 in truth[mintruth][Trace.INFECTED]:
            for node2 in sol[minsol][Trace.INFECTED]:
                if node1 == node2:
                   matchG.add_edge((node1,"truth"),(node2,"sol"),weight=1.0) 
                   continue  
                wscore = 0.0
                dists = set([0.0])
                for first,second in [(node1,node2),(node2,node1)]:
                    try:
                       dist = nx.shortest_path_length(G,source=first,target=second)
                    except nx.NetworkXNoPath:
                       pass
                    else:  
                       dists.add(1.0-(dist/float(diam))) 
                matchG.add_edge((node1,"truth"),(node2,"sol"),weight=max(dists))
        for node1,node2 in matchG.edges():
            assert matchG[node1][node2]["weight"] <= 1.0 and matchG[node1][node2]["weight"] >= 0
        mate = nx.max_weight_matching(matchG, maxcardinality=True)
        for first in mate.keys():
            second = mate[first]
            assert mate[second] == first 
        assert len(mate.keys())/2 == min(len(truth[mintruth][Trace.INFECTED]),len(sol[minsol][Trace.INFECTED]))
        assert len(set(matchG.nodes()).difference(set(mate.keys()))) == abs(len(truth[mintruth][Trace.INFECTED])-len(sol[minsol][Trace.INFECTED]))
        mscore = sum([matchG[key][mate[key]]["weight"] for key in mate.keys()])/2.0
        for remnode in set(matchG.nodes()).difference(set(mate.keys())):
            maxsim = 0.0
            for neighnode in matchG.neighbors(remnode):
                if matchG[remnode][neighnode]["weight"] >= maxsim:
                   maxsim = matchG[remnode][neighnode]["weight"]
            mscore += maxsim
        return mscore/float(max(len(truth[mintruth][Trace.INFECTED]),len(sol[minsol][Trace.INFECTED])))
    
    @staticmethod
    def getPredictionScore(truth,sol):
        """initial spreaders prediction score
        Args:
          truth: truth dict
          sol: our solution dict
        Returns:
          scores: prediction scores 
        """
        state = Trace.INFECTED
        scores = {"sen":0.0,"fpr":0.0,"recall":0.0,"precision":0.0,"acc":0.0,"spec":0.0,"jaccardsim":0.0,"jaccarddis":0.0,"F1":0.0,"F2":0.0,"F01":0.0}
        mintruth = min([time for time in truth.keys() if truth[time].has_key(state)])
        minsol = min([time for time in sol.keys() if sol[time].has_key(state)])
        fn = len(set(truth[mintruth][state]).difference(sol[minsol][state]))
        tp = len(set(truth[mintruth][state]).intersection(sol[minsol][state]))
        fp = len(set(sol[minsol][state]).difference(truth[mintruth][state]))
        union = len(set(sol[minsol][state]).union(truth[mintruth][state]))
        posinf = set()
        for time in truth.keys():
            posinf |= set(truth[time][Trace.INFECTED])
            if truth[time].has_key(Trace.RECOVERED):
               posinf |= set(truth[time][Trace.RECOVERED])
        tn = len(posinf) - tp - fp - fn
        
        scores["sen"] += float(tp) / (tp+fn)
        scores["fpr"] += float(fp) / (fp+tn)
        scores["recall"] += scores["sen"]
        if (tp+fp) != 0:
           scores["precision"] += float(tp) / (tp+fp)
        scores["acc"] += float(tp+tn) / (tp+tn+fn+fp)
        scores["spec"] += 1.0 - scores["fpr"] #true negative rate
        scores["jaccardsim"] = float(tp) / union
        scores["jaccarddis"] = 1.0 - scores["jaccardsim"]
        if scores["precision"] + scores["recall"] != 0.0:
           scores["F1"] += float(2.0*scores["precision"]*scores["recall"]) / (scores["precision"] + scores["recall"])
           scores["F2"] += float(5.0*scores["precision"]*scores["recall"]) / ((4.0*scores["precision"]) + scores["recall"])
           scores["F01"] += float(1.01*scores["precision"]*scores["recall"]) / ((0.01*scores["precision"]) + scores["recall"])
        return scores

    @staticmethod
    def getFootruleScore(ltruth,lest,states=[Trace.INFECTED]):
        """returns partial Spearman's footrule score
        Args:
           ltruth,lest: two dicts
           states: states 
        Returns:
          score: partial footrule score
        """
        PARTIAL = True
        score,INFTIME = 0.0, max(ltruth.keys()) + 1
        truthdif,estdif = max(ltruth.keys()) - min(ltruth.keys()), max(lest.keys()) - min(lest.keys())
        assert truthdif > 0 or estdif > 0
        for state in states:
            nodetimes1 = {node: time for time in ltruth.keys() for node in ltruth[time][state]}
            nodetimes2 = {node: time for time in lest.keys() for node in lest[time][state]}
            allnodes = set(itertools.chain(*[ltruth[time][state] for time in ltruth.keys() if ltruth[time].has_key(state)]))
            for node in allnodes:
                if not nodetimes2.has_key(node):
                   nodetimes2[node]= INFTIME   
            tscore = sum([abs(nodetimes1[node] - nodetimes2[node]) for node in allnodes])
            if PARTIAL:  
               DENOM = float(len(allnodes)) * max(truthdif,estdif)
            else: 
               DENOM = math.ceil(0.5 * len(allnodes)**2)
            score += tscore / DENOM
        return score / len(states)

    @staticmethod
    def getKenTauCorrelationScore(ltruth,lest,states=[Trace.INFECTED]):
        """returns partial Kendall Tau correlation score(kendall b score)
        Args:
          ltruth: truth dict
          lest: estimated dict
          states: states over score will be estimated
        Returns:
          score: partial kendall tau score
        """
        #minest = min(lest.keys())
        #mintruth = min(ltruth.keys())
        #assert mintruth <= 0
        #assert (mintruth < 0 and minest <= 0) or  (mintruth == 0 and minest>=0)
            
        #mintime = min(ltruth.keys())
        #if mintime < 0:
        #   lest2,ltruth2 = deepcopy(lest), deepcopy(ltruth)
        #   lest,ltruth = {}, {} 
        #   for time in ltruth2.keys():
        #       ltruth[time-mintime] = deepcopy(ltruth2[time])
        #   for time in lest2.keys():
        #       lest[time-mintime] = deepcopy(lest2[time]) 
        score,MAXTIME = 0.0, max(ltruth.keys()) + 1
        #assert sorted(lest.keys())[0] >= 0 and sorted(ltruth.keys())[0] >= 0 
        for state in states:
            allnodes = list(set(node for time in ltruth.keys() for node in ltruth[time][state]))
            node2pos = {allnodes[index]:index for index in xrange(len(allnodes))}
            truelist, estlist = [MAXTIME]*len(allnodes), [MAXTIME]*len(allnodes)
            for time in ltruth.keys():
                for node in ltruth[time][state]:
                    truelist[node2pos[node]] = time
            for time in lest.keys():
                for node in set(lest[time][state]).intersection(set(allnodes)):
                    if node2pos.has_key(node):
                       estlist[node2pos[node]] = time
                    
            n0 = (len(allnodes)*(len(allnodes)-1))/2 
            con,dis,n1,n2 = 0,0,0,0
            for index1 in xrange(len(allnodes)):
                pos1 = node2pos[allnodes[index1]]
                for index2 in xrange(index1+1,len(allnodes)):
                    pos2 = node2pos[allnodes[index2]]
                    if (estlist[pos1] - estlist[pos2]) * (truelist[pos1] - truelist[pos2]) > 0:
                       con += 1 
                    if (estlist[pos1] - estlist[pos2]) * (truelist[pos1] - truelist[pos2]) < 0:
                       dis += 1 
                    if estlist[pos1] == estlist[pos2]: 
                       n1 += 1
                    if truelist[pos1] == truelist[pos2]: 
                       n2 += 1  
            if n0 == n1 or n0 == n2:
               denom = math.sqrt((n0-n1+1)*(n0-n2+1))
            else:
               denom = math.sqrt((n0-n1)*(n0-n2))             
            kenbscore = float(con-dis)/denom
            score += kenbscore  #print sp.stats.kendalltau(truelist, estlist)[0]
        return score / len(states)
    
      
    @staticmethod
    def getKenTauScore(ltruth,lest,p=0.5,states=[Trace.INFECTED]):
        """returns partial Kendall Tau score
        Args:
          ltruth: truth dict
          lest: estimated dict
          p: p parameter for partial ranking
          states: states over score will be estimated
        Returns:
          score: partial kendall tau score
        """
        return DiscreteHistoryScore.getInternalKendallScore(ltruth,lest,p,states)
    
    @staticmethod
    def getHausdorffScore(ltruth,lest,states=[Trace.INFECTED]):
        """returns Hausdorff score
        Args:
           ltruth:
           lest:
           states:
        Returns:
           score:  
        """
        return DiscreteHistoryScore.getInternalKendallScore(ltruth,lest,"Hausdorff",states)
        
    @staticmethod
    def getInternalKendallScore(ltruth,lest,param,states):
        """returns partial Kendall Tau score or Hausdorff score
        Args:
          ltruth: truth dict
          lest: estimated dict
          param: p parameter for partial ranking for Kendall Tau or Hausdorff
          states: states over score will be estimated
        Returns:
          score: partial kendall tau score
        """
        sortedtimes1, sortedtimes2 = sorted(ltruth.keys()), sorted(lest.keys())
        score, INFTIME = 0.0, 100000
        for state in states:
            truthtimes,esttimes = {}, {}
            for time in sortedtimes1:
                if ltruth[time].has_key(state):
                   for node in ltruth[time][state]: 
                       truthtimes[node] = time
            for time in sortedtimes2:
                if lest[time].has_key(state):
                   for node in lest[time][state]:  
                       esttimes[node] = time       
            allnodes = truthtimes.keys()
            for node in allnodes:
                if not esttimes.has_key(node):
                   esttimes[node] = INFTIME   
            U, S, T = 0.0, 0.0, 0.0
            for index1 in xrange(len(allnodes)):
                node1 = allnodes[index1]
                for index2 in xrange(index1+1,len(allnodes)):
                    node2 = allnodes[index2]
                    if (truthtimes[node1] - truthtimes[node2]) == 0 and (esttimes[node1] - esttimes[node2]) != 0:
                       S += 1
                    elif (esttimes[node1] - esttimes[node2]) == 0 and (truthtimes[node1] - truthtimes[node2]) != 0:
                       T += 1
                    elif (truthtimes[node1] - truthtimes[node2]) * (esttimes[node1] - esttimes[node2]) < 0:
                       U += 1
            if param == "Hausdorff":           
               curscore =  U + max(S,T)   
            elif type(param) == float:
               curscore = U + (param * (S + T))
            if len(allnodes) > 1:    
               curscore /= float(0.5* len(allnodes) * (len(allnodes) - 1))
            else:
               curscore = 1.0 
            score += curscore
        return score / len(states)
