import sys
sys.path.append("./lib")
import os
import time
import HistoryInferRunner
from Trace import Trace

def assignDefaultHistoryParameters():
    """assigns default parameter values
    Returns:
       returns default value array
    """
    vararr = {}
    vararr["tracefile"] = "tracefile.trace"
    vararr["snapshotfile"] = "trial.snapshot"
    vararr["resultfile"] = "histresults.hist"
    vararr["runfolder"] = "temprun"
    vararr["smodel"] = "si"
    vararr["evol"] = "static"
    vararr["algo"] = "GreedySubCoverSingle"
    vararr["algopar"] = {}
    vararr["noisetype"] = None
    vararr["infermode"] = "Spreader"
    vararr["dist"] = "cont"
    vararr["graphfile"] = "trial.edgelist" 
    vararr["printscore"] = "KenTau"
    vararr["scoretimes"] = None
    vararr["scorestates"] = "{0}".format(Trace.INFECTED)
    return vararr

def readHistoryConfigFile(filename):
    """reads CEFER configuration file and returns parameters
    Args:
       filename: configuration filename
    Returns:
       returns parameters from configuration file
    """
    vararr = assignDefaultHistoryParameters()
    with open(filename,"r") as file:
        for line in file:
            varname,value = line.rstrip().split(":")
            value = value.lstrip()
            if varname in ["iter","objmethod","roundmethod","approx","method"]:
               vararr["algopar"][varname] = value
            elif varname in ["ensemble"]:
               vararr["algopar"][varname] = value == "True"
            elif varname in ["enscount"]:
               vararr["algopar"][varname] = int(value)      
            elif varname in ["noise"]:
               vararr[varname] = float(value)
            else:
               vararr[varname] = value
    assert vararr["smodel"] in ["si","sir","seir","sis"] 
    return vararr


def main():
   """runs arbitrary HISTORY inference by the parameters provided in config file and outputs the inferred graph
   Args:.
      configfile: configuration filename
   """
   if len(sys.argv) != 2:
      print "Configuration file must be given\n Look at example.config for format details!!"
      exit(1)
   vararr = readHistoryConfigFile(sys.argv[1])
   t1 = time.time()
   HistoryInferRunner.runHistory(vararr)
   t2 = time.time()
   print "Total history inference done in {0} seconds".format(t2 - t1)
   
if  __name__ == '__main__':
    main()

    
