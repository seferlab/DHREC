import os
import sys
sys.path.append("./lib")
import myutilities as myutil
import math
import time
from copy import deepcopy
sys.path.append("./histcodes")
import DiscreteBound
import DiscreteGreedy
from Trace import Trace
from InputOutput import InputOutput
from DiscreteHistoryScore import DiscreteHistoryScore as DisScore
import QsapArbitraryExternal
import QsapLogconcaveExternal
import IterPCDSVC
import IterMinCut
import Independent
import greedyForward

def getScore(resultfile,score,tracefile,scoretimes,inter,evol,infermode,scorestates,G,smodel,dist):
    """ estimates score of resultfile
      Args:
         resultfile:
         score:
         tracefile:
         scoretimes:
         inter:
         evol:
         infermode:
         scorestates:
         G:
         smodel:
         dist:
      Returns:
         score: score
    """
    if infermode == "Spreader":
       assert score not in ["KenTau","Footrule","Hausdorff","KenTauCorrelation"]
    elif infermode == "History":
       assert score not in ["Graph","F1","F01","jaccardsim","jaccarddis"]
    trace = InputOutput.readPlainTrace(tracefile,dist)
    truth = DisScore.trace2History(trace,evol,smodel,scoretimes) #exact state transition times not snapshots
    mintime = min(truth.keys())
    if inter == "bound":
       if mintime < 0:
          truth2 = deepcopy(truth)
          truth = {}
          for time in truth2.keys():
              truth[time-mintime] = deepcopy(truth2[time])
              
    if infermode == "Spreader":
       mintime,sol = InputOutput.readHistoryResult(resultfile,infermode)
       sol = {mintime: sol}
    elif infermode == "History":
       sol = InputOutput.readHistoryResult(resultfile,infermode)
    for scorestate in scorestates:
        for time in truth.keys():
            if not truth[time].has_key(scorestate):
               truth[time][scorestate] = set() 
        for time in sol.keys():
            if not sol[time].has_key(scorestate):
               sol[time][scorestate] = set()  

    splitted = score.split(",")
    if len(splitted) != 1:
       score = splitted[1]
       methodname = "getPartGraphScore"
       return getattr(DisScore,methodname)(truth,sol,G,10,score,smodel,evol,scoretimes,inter)
    if score == "KenTau":  
       methodname = "get{0}Score".format(score)
       return getattr(DisScore,methodname)(truth,sol,0.5,scorestates)
    elif score in ["Footrule","Hausdorff","KenTauCorrelation"]: 
       methodname = "get{0}Score".format(score)
       return getattr(DisScore,methodname)(truth,sol,scorestates)  
    elif score in ["sen","fpr","recall","precision","acc","spec","F1","F01","F01","jaccardsim","jaccarddis"]:
       methodname = "getPredictionScore"
       return getattr(DisScore,methodname)(truth,sol)[score]
    
def runChecks(G,algo,dist,inter,smodel):
    """checks for parameter prerun
    Args:
       G: 
       algo: 
       dist:
       inter:
       smodel:
    Returns:
       bool: True/False
    """
    if dist == "dis" and inter == "bound":
       if algo not in ["second-logconcave-external","second-arbitrary-external","second-internal","Independent","MatroidSub"]:
          print "Algo {0} is unknown!".format(algo) 
       assert algo in ["second-logconcave-external","second-arbitrary-external","second-internal","Independent","MatroidSub"]
    elif dist == "dis" and inter == "None":
       if algo not in ["Pcvc","Pcdsvc","MinCut","Independent","GreedySubCoverSingle","GreedyForward"]:
          print "Algo {0} is unknown!".format(algo)   
       assert algo in ["Pcvc","Pcdsvc","MinCut","Independent","GreedySubCoverSingle","GreedyForward"] 
       if algo in ["Pcvc","Pcdsvc"]:
          assert smodel in ["si","sir"]
       elif algo == "MinCut":
          assert smodel == "seir"
    if dist == "dis":
       return True 
       for edge in G.edges(data=True):
           if edge[2][Trace.SPROB] != 1.0:
              return False 
       return True    

def runHistory(vararr):
    """ infers history
    Args:
       vararr: arguments in dict
    """
    for var in vararr.keys():
        if type(vararr[var]) == type(""):
           exec('{0}="{1}"'.format(var,vararr[var]))
        else:
           exec('{0}={1}'.format(var,vararr[var]))   
    scorestates = scorestates.split(" ")
    curstates,obstimes = InputOutput.readSnapshot(snapshotfile,dist)
    G = InputOutput.readGraphAndParams(graphfile)
    assert min(obstimes) == 0 and len(curstates) == len(obstimes) and runChecks(G,algo,dist,inter,smodel) 
    algoparstr = "\n".join(["{0} {1}".format(key,algopar[key]) for key in algopar.keys()])
    print "Inferring History by the algorithm {0} with parameters {1}".format(algo,algoparstr)
    print "on {0} {1} case from {2} observed points".format(dist,inter,len(obstimes))
    print "Using the observed points {0}, 0 means the earliest time point".format(obstimes)
    print "Using snapshots from file {0}".format(snapshotfile)
    stime = time.time()    
    if dist == "dis" and inter == "None":
       if algo in ["Pcvc","Pcdsvc"]:
          IterPCDSVC.runner(algo,algopar,curstates,G,smodel,infermode,resultfile,runfolder,obstimes,True)
       elif algo == "MinCut":
          IterMinCut.runner(algopar,curstates,G,smodel,infermode,resultfile,runfolder,True) 
       elif algo == "Independent":
          Independent.runner(inter,G,smodel,curstates,obstimes,resultfile,True) 
       elif algo == "GreedyForward":
          greedyForward.runner(inter,G,smodel,curstates,obstimes,resultfile,True)       
       elif algo == "GreedySubCoverSingle":
          DiscreteGreedy.runner(algo,algopar,infermode,curstates,G,smodel,resultfile,runfolder,obstimes,True)    
    elif dist == "cont" and inter == "None":
       ContinuousGreedy.runner(algo,algopar,infermode,curstates,G,smodel,resultfile,runfolder)
    endtime = time.time()
    print "Inference took {0} seconds".format(endtime-stime)    
    if printscore != "None":
       if inter == "None":
          assert scoretimes != "None"
       elif inter == "bound":
          scoretimes = obstimes
       score = getScore(resultfile,printscore,tracefile,scoretimes,inter,evol,infermode,scorestates,G,smodel,dist)
       print "{0} Prediction score {1}: {2} ".format(infermode,printscore,score)
