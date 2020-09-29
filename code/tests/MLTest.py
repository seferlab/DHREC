import random
import unittest
import sys
import networkx as nx
import autohistoryinfer
import autohistorytracegen

class HistoryInferTest(unittest.TestCase):
    """ Tests for QSAP History Reconstruction
    """
    
    def setUp(self,testname,algo,algopar,G,tracepath,obstimes,noise,startcount,timefrac,smodel,infermode,noisetype,resultpath,result,sol):
        """sets up test
        Args:
           testname: 
           algo:
           algopar:
           G: graph with dist attributes
           tracepath: trace path
           obstimes: observed times
           noise: noise 
           startcount: start node count
           timefrac: timefrac
           smodel: spreading model
           infermode: infermode
           noisetype:
           resultpath:
           result:
           sol:
       """
       super(HistoryInferTest, self).__init__(testname)
       self.G = G
       self.algo = algo
       self.algopar = algopar
       self.tracepath = tracepath
       self.obstimes = obstimes
       self.noise = noise
       self.startcount = startcount
       self.timefrac =  timefrac
       self.smodel = smodel
       self.infermode = infermode
       self.noisetype = noisetype
       self.resultpath = resultpath
       self.result = result 
       self.curstates = curstates
       if self.infermode == "History":
          self.itimes = {}
          self.rtimes = {} 
          self.stimes = {} 
       self.sol = sol   
        
    def checkVariableValidity(modcurstates,varname,smodel):
        """check whether variable is valid on given modcurstates
        Args:
           modcurstates:
           varname:
           smodel:
        """
        sortedtimes = sorted(modcurstates.keys())
        nodetimes = var2Times(varname)
        if smodel == "si":
           for time in sortedtimes:
               if time < nodetimes[Trace.INFECTED]:
                  assert node in modcurstates[time][Trace.SUSCEPTIBLE]
              else:
                  assert node in modcurstates[time][Trace.INFECTED] 
        elif smodel == "sir":
           if nodetimes.has_key(Trace.INFECTED) and nodetimes.has_key(Trace.RECOVERED):
              for time in sortedtimes:
                  if time < nodetimes[Trace.INFECTED]:
                     assert node in modcurstates[time][Trace.SUSCEPTIBLE]
                  elif time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                     assert node in modcurstates[time][Trace.INFECTED]
              elif time >= nodetimes[Trace.RECOVERED]:
                 assert node in modcurstates[time][Trace.RECOVERED]    
           elif nodetimes.has_key(Trace.INFECTED):
          for time in sortedtimes:
              if time < nodetimes[Trace.INFECTED]:
                 assert node in modcurstates[time][Trace.SUSCEPTIBLE]
              else:
                 assert node in modcurstates[time][Trace.INFECTED]
              #assert node not in modcurstates[time][Trace.RECOVERED]         
    elif smodel == "seir":
       if nodetimes.has_key(Trace.EXPOSED) and nodetimes.has_key(Trace.INFECTED) and nodetimes.has_key(Trace.RECOVERED):
          assert nodetimes[Trace.EXPOSED] > 0 
          for time in sortedtimes:
              if time < nodetimes[Trace.EXPOSED]:
                 assert node in modcurstates[time][Trace.SUSCEPTIBLE]
              elif time >= nodetimes[Trace.EXPOSED] and time < nodetimes[Trace.INFECTED]:
                 assert node in modcurstates[time][Trace.EXPOSED]
              elif time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                 assert node in modcurstates[time][Trace.INFECTED]   
              elif time >= nodetimes[Trace.RECOVERED]:
                 assert node in modcurstates[time][Trace.RECOVERED]   
       elif nodetimes.has_key(Trace.INFECTED) and nodetimes.has_key(Trace.RECOVERED):
          assert nodetimes[Trace.INFECTED] == 0  
          for time in sortedtimes:
              elif time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                 assert node in modcurstates[time][Trace.INFECTED]
              elif time >= nodetimes[Trace.RECOVERED]:
                 assert node in modcurstates[time][Trace.RECOVERED]       
       elif nodetimes.has_key(Trace.EXPOSED) and nodetimes.has_key(Trace.INFECTED):
          assert nodetimes[Trace.EXPOSED] > 0 
          for time in sortedtimes:
              elif time >= nodetimes[Trace.INFECTED] and time < nodetimes[Trace.RECOVERED]:
                 assert node in modcurstates[time][Trace.INFECTED]
              elif time >= nodetimes[Trace.RECOVERED]:
                 assert node in modcurstates[time][Trace.RECOVERED] 
       elif nodetimes.has_key(Trace.INFECTED):
          assert nodetimes[Trace.INFECTED] == 0 
       elif nodetimes.has_key(Trace.EXPOSED):
          assert nodetimes[Trace.EXPOSED] > 0 
    return True

          
    def testHistoryValidity(self):
        """tests validity of history
        Args:
        Returns:
        """
        if self.infermode == "History":
           sortedtimes = sorted(self.result.keys())
           for index1 in xrange(1,len(sortedtimes)):
               time1 = sortedtimes[index1]
               for index2 in xrange(index1):
                   time2 = sortedtimes[index2]
                   for node1 in self.result[time1][Trace.INFECTED]: 
                       validnodes = set(node for node in self.result[time2][Trace.INFECTED] if not rtimes.has_key(node) or (rtimes.has_key(node) and rtimes[node] >= time1))
                       self.AssertGreater(len(set(G.predecessors(node1)).intersection(validnodes)),0)
        elif self.infermode == "Spreader":
           mintime = min(self.result.keys())
           inodeset = set(self.result[mintime][Trace.INFECTED])
           self.assertEqual(len(inodeset.difference(set(self.curstates[0][Trace.INFECTED] + self.curstates[0][Trace.RECOVERED]))),0) 

           
def getTests(smodel,noise,prob,G,distinfo,traces):
    """generates tests for discrete none case
    Args:
       smodel:
       noise:
       prob:
       G:
       distinfo:
       traces:
    Returns:
    """
    suiteFew = unittest.TestSuite()
    suiteFew.addTest(DiscreteNoneTest("testHistoryValidity",traces,smodel,distinfo,G,noise,samplerate,prob))
    return suiteFew    
    
def getAlgos():
    """
    """
    return

def main():
    """ main function for testing traces
    """
    graphpref = "graphs"
    newtracepref = "testtracesnapshots"
    tracepref = "testtraces"
    configpref = "testhistconfig"
    runpref = "testrun"
    resultpref = "testresult"
    resultfolderpref = "."
    runfolderpref = "."
    tracefolderpref = ".."
    newtracefolderpref = "."
    graphfolderpref = ".."
    configfolderpref = "."
    pbsfolder = "pbsfolder"
    realdatas = ["real","syn"]
    evols = ["static"]
    infermodes = ["History","Spreader"]
    smodels = ["si","sir","seir"]
    probs = ["dis","cont"]
    inters = ["bound",None] 
    startcounts = [1,3,5,6,8,10,12] #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    timefracs = [0.51,0.8] #percentage of diffusion maxtime(peak time gibi) [0.51,0.6,0.7,0.8,0.9]
    timecounts = [1,2,3,4,5] 
    noises = [0.0,0.2,0.4,0.6] #[0.0,0.2,0.4,0.6,0.8,0.9]
    samplerates = [0]
    noisetype = "StateChange" #Normal noise
    filesamplecount = 10 #for syn data
    fraccons = {Trace.INFECTED: {"min":0.001, "max":1.0}}
    printscore = "Hausdorff"  #"Hausdorff" "Graph,Prediction,F1" #??
    tracecount = 2
    sismaxtime = 100
    
    parts = list(itertools.product(noises,realdatas,evols,timefracs,timecounts,smodels,startcounts,infermodes,inters,samplerates,probs))
    random.shuffle(parts)
    for (noise,realdata,evol,timefrac,timecount,smodel,startcount,infermode,inter,samplerate,prob) in parts:
        if inter == None and timecount != 1:
           continue
        
        (graphfolder,tracefolder,configfolder) = autohistorytracegen.genMainFolders(graphpref,graphfolderpref,tracepref,tracefolderpref,configpref,configfolderpref,pbsfolder,realdata,smodel,evol,prob)
        path2info = autohistorytracegen.returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount,prob,smodel)
        for path in path2info.items():
            path = item[0]
            Gname,filetracefolder = item[1]

            #Trace gen
            extensionstr = "-".join(path.replace("../","").split("/")[1:])
            noisesamplestr = "{0}-{1}-{2}".format(noisetype,noise,samplerate)
            tracefolder = "{0}/{1}/{2}".format(filetracefolder,noisesamplestr,startcount)

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
            vararr["dists"] = {}
            vararr["dist"] = prob
            vararr["modelparams"] = {}
            
            if smodel == "sis":
               vararr["sismaxtime"] = sismaxtime
             
            configpath = "{0}/{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}-{9}.config".format(configfolder,extensionstr,realdata,evol,smodel,prob,tracecount,samplerate,noise,startcount)
            genTraceConfigFile(configpath,vararr)
            code = "python tracegen.py {0}".format(configpath)
            os.system(code)
            continue
            if random.random() <= 0.5:
               POOL = "pool1"
            else:
               POOL = "pool2"  
            pbsfilename="{0}-{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}.pbs".format(realdata,evol,smodel,prob,extensionstr,tracecount,samplerate,noise,startcount)
            Pbs.submitPbs(code,pbsfolder,pbsfilename,POOL)
            
            #Inference Part
            graphtraceinput = [graphfolderpref,tracepref,tracefolderpref,newtracepref,newtracefolderpref]
            runresultinput = [runpref,runfolderpref,resultpref,resultfolderpref]
            configinput = [configpref,configfolderpref,pbsfolder]
            (graphfolder,tracefolder,newtracefolder,configfolder,resultfolder,runfolder) = genMainFolders(graphtraceinput,runresultinput,configinput)
               
            #graphinput = [graphfolder,realdata,evol,prob,smodel,filesamplecount,inter]
            #traceinput = [tracefolder,fraccons,samplenoisestr,startcount]
            #path2info = returnPathInfo(graphinput,traceinput)
            pathkeys = autohistoryinfer.path2info.keys()
            random.shuffle(pathkeys)
            path = pathkeys[0]
            algoblocks = getAlgoBlocks(prob,inter,infermode,smodel)
            random.shuffle(algoblocks)
            algoblock = algoblocks[0]
            Gname,tracefile = path2info[path]
            tracestr = "-".join(tracefile.split("/")[-3:])
            algo,algoparams = algoblock
            
            algofields = [algo,algoparams,infermode,inter,runfolder]
            tracefields = [noise,noisetype,timefrac,timecount,prob,smodel,evol,realdata,tracefile,tracestr,newtracefolder]
            otherfields = [complete,completeupto,path,printscore]
            autohistoryinfer.runMyAlgos(algofields,tracefields,otherfields)
            resultfilename = ""
            
            noisesamplestr = "{0}-{1}-{2}".format(noisetype,noise,samplerate)
            traceoutfolder = "{0}/{1}/{2}".format(filetracefolder,noisesamplestr,startcount)
            (traces,allnodes) = InputOutput.readTraces(traceoutfolder,tracecount,smodel)
            [os.system("rm -rf {0}".format(folder)) for folder in [traceoutfolder,configpath]]
            suiteFew = getTests(smodel,noise,prob,G,spreadparam,traces,sol)
            unittest.TextTestRunner(verbosity = 2).run(suiteFew)
                                  
       

if __name__ == '__main__':
    main()
       
