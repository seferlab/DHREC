import networkx as nx
import sys
sys.path.append("./lib")
import os
import gzip
import cPickle
import random
import myutilities as myutil
import itertools
from InputOutput import InputOutput
from Trace import Trace
from Pbs import Pbs
import HistoryInferRunner
import HistOtherAlgos

def genOtherConfigFile(configfile,graphfile,algo,algoparams,snapshotfile,resultfile,smodel,trunfolder,inter):
    """generates config file for other algos
    Args:
       configfile:
       config parameters:
    """
    with gzip.open(configfile,"wb") as file:
        cPickle.dump(graphfile,file)
        cPickle.dump(algo,file)
        cPickle.dump(algoparams,file)
        cPickle.dump(snapshotfile,file)
        cPickle.dump(resultfile,file)
        cPickle.dump(smodel,file)
        cPickle.dump(trunfolder,file)
        cPickle.dump(inter,file)
 
def genHistoryConfigFile(configfile,vararr):
    """generates history config file
    Args:
       configfile:
       vararr:
    """
    with open(configfile,"w") as file:
        for var in vararr.keys():
            if var == "dists":
               for key in vararr[var].keys():
                   if key in [Trace.S2I, Trace.I2R, Trace.E2I, Trace.I2S, Trace.S2E]:
                      paramstr = " ".join([str(item) for item in vararr[var][key][1]])
                      file.write("{0}: {1} {2}\n".format(key,vararr[var][key][0],paramstr))
                   elif key == Trace.SPROB:
                      file.write("{0}: {1}\n".format(key,vararr[var][key])) 
            else:
               file.write("{0}: {1}\n".format(var,vararr[var]))

def getAlgoBlocks(prob,inter,infermode,smodel):
    if prob == "dis" and inter == "bound":
       if infermode == "Spreader":
          algoblocks = [("MatroidSub",{"method":"search","ensemble":False}),("MatroidSub",{"method":"search","ensemble":True,"enscount":2}),("MatroidSub",{"method":"search","ensemble":True,"enscount":5})]
          algoblocks = [("second-logconcave-internal",{"ensemble":False,"approx":"reliability"}),("second-arbitrary-internal",{"ensemble":False,"method":"greedy"})]
          algoblocks = [("second-arbitrary-internal",{"ensemble":False,"method":"pipage"})]
          algoblocks.extend(algoblocks2)
          otheralgos = [("RumorCentrality",{})]
          algoblocks.extend(otheralgos)
       elif infermode == "History":
          algoblocks = [("MatroidSub",{"method":"search","ensemble":False})]
          #algoblocks = [("second-logconcave",{"ensemble":False,"approx":"reliability"}),("second-arbitrary",{"ensemble":False,"method":"greedy"})]
          #algoblocks = [("second-arbitrary",{"ensemble":False,"method":"pipage"})]
          #algoblocks = [("MatroidSub",{"method":"search","ensemble":False})]
          #algoblocks = [("Qsapmin",{"ensemble":False,"approx":"reliability"})]
          #algoblocks = [("Qsapmax",{"ensemble":False,"method":"pipage"})]
          #algoblocks.extend(algoblocks2)
          #indeblocks = [("Independent",{"ensemble":False})] 
          #algoblocks.extend(indeblocks)  
    elif prob == "dis" and inter == None:
       if infermode == "Spreader":
          algoblocks = [("GreedySubCoverSingle",{"iter":None,"objmethod":"log","ensemble":False}),("GreedySubCoverSingle",{"iter":None, "objmethod":"log","ensemble":True,"enscount":2})]
          #,("GreedySubCoverSingle",{"iter":None, "objmethod":"log","ensemble":True,"enscount":3})
          #algoblocks = [("FracCover",{"iter":None, "objmethod":"log","ensemble":False,"roundmethod":"random"})]
          #algoblocks = []

          if smodel in ["si","sir"]:
             appralgos = [("Pcdsvc",{"ensemble":False}),("Pcdsvc",{"ensemble":True,"enscount":2})]
             appralgos2 = [("Pcvc",{"ensemble":False}),("Pcvc",{"ensemble":True,"enscount":2})]
             algoblocks.extend(appralgos)
             algoblocks.extend(appralgos2)
          elif smodel == "seir":
             appralgos = [("MinCut",{"ensemble":False}),("MinCut",{"ensemble":True,"enscount":2})]       
             algoblocks.extend(appralgos)
             
          #algoblocks = []
          #otheralgos = [("NetSleuth",{}),("RumorCentrality",{}),("KEffectors",{})]
          #otheralgos = [("NetSleuth",{})]
          otheralgos = [("RumorCentrality",{})]
          algoblocks.extend(otheralgos)
       elif infermode == "History":
          algoblocks = [("GreedySubCoverSingle",{"iter":None, "objmethod":"log","ensemble":False})]
          
          if smodel in ["si","sir"]:
             appralgos = [("Pcdsvc",{"ensemble":False}),("Pcvc",{"ensemble":False})]
          elif smodel == "seir":
             appralgos = [("MinCut",{"ensemble":False})]
          algoblocks.extend(appralgos)    
          
          indeblocks = [("GreedyForward",{"ensemble":False})] 
          algoblocks.extend(indeblocks)
    elif prob == "cont" and inter == "bound":
       algoblocks = [("greedy",{})]  
    elif prob == "cont" and inter == None:
       algoblocks = [("greedy",{})]  
    return algoblocks

def returnAlgoStr(algo,algoparams):
    return "-".join([algo] + [str(item) for item in algoparams.values()])

def genMainFolders(graphtraceinput,runresultinput,configinput,sideinput):
    (graphfolderpref,tracepref,tracefolderpref,newtracepref,newtracefolderpref) = graphtraceinput
    (runpref,runfolderpref,resultpref,resultfolderpref) = runresultinput
    (configpref,configfolderpref,pbsfolder) = configinput    
    (smodel,prob,evol,realdata,inter) = sideinput 
    graphfolder = "{0}/{1}_{2}_{3}".format(graphfolderpref,realdata,evol,"graphs")
    tracefolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(tracefolderpref,tracepref,realdata,evol,smodel,"edge",prob)
    newtracefolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(newtracefolderpref,newtracepref,realdata,evol,smodel,prob,inter)
    configfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(configfolderpref,configpref,realdata,evol,smodel,prob,inter)
    resultfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(resultfolderpref,resultpref,realdata,evol,smodel,prob,inter)
    runfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(runfolderpref,runpref,realdata,evol,smodel,prob,inter)
    [os.makedirs(folder) for folder in [pbsfolder,configfolder,resultfolder,runfolder] if not os.path.exists(folder)]
    return graphfolder,tracefolder,newtracefolder,configfolder,resultfolder,runfolder


def returnPathInfo(graphinput,traceinput):
    """returns path info
    Args:
       graphinput:
       traceinput:
    Returns:
       path2info:   
    """
    (graphfolder,realdata,evol,prob,smodel,filesamplecount,inter) = graphinput
    (tracefolder,fraccons,samplenoisestr,startcount) = traceinput
    path2info={}
    if realdata == "real" and evol == "static":
       for filename in myutil.listfiles(graphfolder):
           if filename.split("-")[0:2] != [prob,smodel]:
              continue
           filepath ="{0}/{1}".format(graphfolder,filename)
           if inter == None and (filename.find("weibull") != -1 or filename.find("rayleigh") != -1 or filename.find("powerlaw") != -1):
              continue
           #if filename.find("sn") != -1:
           #   if filename.find("sn1") == -1: 
           #      continue
           #if filename.find("grid") != -1:
           #   continue     
           filetracefolder = "{0}/{1}/{2}/{3}".format(tracefolder,filename,samplenoisestr,startcount)
           assert os.path.exists(filetracefolder)
           G = InputOutput.readGraphAndParams(filepath)
           minicount = int(round(G.number_of_nodes() * fraccons[Trace.INFECTED]["min"]))
           maxicount = int(round(G.number_of_nodes() * fraccons[Trace.INFECTED]["max"]))
           sentfraccons = {Trace.INFECTED: {"min":minicount, "max":maxicount}}
           while True:
               filem = Trace.getRandomTraceFile(filetracefolder,sentfraccons)   
               if filem not in ["0.plain","1.plain"]:
                  continue
               tracefile = "{0}/{1}".format(filetracefolder,filem)
               break 
           #tracefile = "{0}/{1}".format(filetracefolder,Trace.getRandomTraceFile(filetracefolder,sentfraccons))
           path2info[filepath] = (filename,tracefile)
    elif realdata == "syn" and evol == "static":
       for filename in myutil.listfiles(graphfolder):
           filepath = "{0}/{1}".format(graphfolder,filename)
           if filename.split("-")[0:2] != [prob,smodel]:
              continue
           if inter == None and (filename.find("weibull") != -1 or filename.find("rayleigh") != -1 or filename.find("powerlaw") != -1):
              continue
           innum = int(filename.split("_")[-1].replace(".edgelist",""))
           if innum > 1:
              continue
           if smodel == "si":
              if filename.find("sprob_1.0_1.0") != -1 and (filename.find("expo_0.2_1.0") != -1 or filename.find("expo_0.5_0.5") != -1 or filename.find("expo_0.1_0.5") != -1):
                 pass
              else:
                 continue
           elif smodel == "sir":
              if filename.find("sprob_1.0_1.0") != -1 and filename.find("s2i_expo_0.2_1.0") != -1 and filename.find("i2r_expo_0.8_1.0") != -1:
                 pass
              else:
                 continue 
           elif smodel == "seir":
              if filename.find("sprob_1.0_1.0") != -1 and filename.find("s2e_expo_0.5_2.0") != -1 and filename.find("e2i_expo_0.5_2.0") != -1 and filename.find("i2r_expo_0.5_2.0") != -1:
                 pass
              else:
                 continue 
           filetracefolder = "{0}/{1}/{2}/{3}".format(tracefolder,filename,samplenoisestr,startcount)
           assert os.path.exists(filetracefolder)
           G = InputOutput.readGraphAndParams(filepath)
           minicount = int(round(G.number_of_nodes() * fraccons[Trace.INFECTED]["min"]))
           maxicount = int(round(G.number_of_nodes() * fraccons[Trace.INFECTED]["max"]))
           sentfraccons = {Trace.INFECTED: {"min":minicount, "max":maxicount}}  
           while True:
               filem = Trace.getRandomTraceFile(filetracefolder,sentfraccons)   
               if filem not in ["0.plain","1.plain"]:
                  continue
               tracefile = "{0}/{1}".format(filetracefolder,filem)
               break 
           #tracefile = "{0}/{1}".format(filetracefolder,Trace.getRandomTraceFile(filetracefolder,sentfraccons))
           path2info[filepath] = (filename,tracefile)
    return path2info

def runOtherAlgos(algofields,tracefields,otherfields):
    """runs other algos
    Args:
       algofields:
       tracefields:
       otherfields:
    Returns:
    """
    return
    algo,algoparams,runfolder = algofields
    noise,noisetype,timefrac,timecount,realdata,tracefile,tracestr,newtracefolder,tracefolder = tracefields
    complete,completeupto,path,graphfolder,configfolder,resultfolder = otherfields
    algostr = returnAlgoStr(algo,algoparams)
    parameterstr = "frac{0}-count{1}".format(timefrac,timecount)
    extensionstr = "-".join(path.replace("./","").split("/")[1:])
    tresultfolder = "{0}/{1}/{2}/{3}/{4}/{5}".format(resultfolder,extensionstr,tracestr,parameterstr,infermode,algostr)
    if not os.path.exists(tresultfolder):
       os.makedirs(tresultfolder)
    indices = set([-1] + [int(myfile.replace(".hist","")) for myfile in myutil.listfiles(tresultfolder)])
    if complete and len(indices) >= completeupto + 1:
       return
    resultfile = "{0}/{1}.hist".format(tresultfolder,max(indices)+1)
    trunfolder = "{0}/{1}/{2}/{3}/{4}/{5}/{6}".format(runfolder,extensionstr,tracestr,parameterstr,infermode,algostr,max(indices)+1)
    if not os.path.exists(trunfolder):
       os.makedirs(trunfolder)
    tnewtracefolder = "{0}/{1}/{2}/{3}/{4}/{5}/{6}".format(newtracefolder,extensionstr,tracestr,parameterstr,infermode,algostr,max(indices)+1)   
    if not os.path.exists(tnewtracefolder):
       os.makedirs(tnewtracefolder)
    snapshotfile = "{0}/infer.snapshot".format(tnewtracefolder)

    G = InputOutput.readGraphAndParams(path)
    assert timecount == 1
    if infermode == "History":
       fracpoints = [(timefrac*index)/timecount for index in xrange(1,timecount+1)]
    elif infermode == "Spreader":
       fracpoints = [timefrac]
       for index in xrange(1,timecount):
           interval = (1.0-timefrac)/(timecount-1)
           fracpoints.append(timefrac+(index*interval))
       assert max(fracpoints) <= 1.00000001
    trace = InputOutput.readPlainTrace(tracefile,prob)
    maxtime = Trace.getMaxTraceTime(trace)
    obstimes = sorted(list(set([int(round(frac*maxtime)) for frac in fracpoints])))
    if 0 in obstimes:
       obstimes.remove(0)
    if len(obstimes) == 0:
       return     
    curstates = [Trace.trace2Snapshot(trace,obstime,smodel,G) for obstime in obstimes]
    InputOutput.writeSnapshot(curstates,infermode,inter,smodel,obstimes,snapshotfile)
        
    configpath = "{0}/config_{1}_{2}_{3}_{4}_{5}-{6}.config".format(configfolder,extensionstr,infermode,algostr,parameterstr,tracestr,max(indices)+1)
    genOtherConfigFile(configpath,path,algo,algoparams,snapshotfile,resultfile,smodel,trunfolder,inter)
    code = "python HistOtherAlgos.py {0}".format(configpath)
    #os.system(code)
    #return
    if random.random() <= 0.5:
       pool = "pool2"
    else:
       pool = "pool2" 
    pbsfilename = "{0}-{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}-{9}-{10}.pbs".format(realdata,evol,smodel,prob,inter,extensionstr,infermode,algostr,parameterstr,tracestr,max(indices)+1)
    Pbs.submitPbs(code,pbsfolder,pbsfilename,pool)

    
def runMyAlgos(algofields,tracefields,otherfields):
    """runs my algos
    Args:
       algofields:
       tracefields:
       otherfields:
    Returns:
    """
    algo,algoparams,infermode,inter,runfolder = algofields
    noise,noisetype,timefrac,timecount,prob,smodel,evol,realdata,tracefile,tracestr,newtracefolder,tracefolder = tracefields
    complete,completeupto,path,printscore,graphfolder,configfolder,resultfolder = otherfields
    algostr = returnAlgoStr(algo,algoparams)
    parameterstr = "frac{0}-count{1}".format(timefrac,timecount)
    extensionstr = "-".join(path.replace("./","").split("/")[1:])
    tresultfolder = "{0}/{1}/{2}/{3}/{4}/{5}".format(resultfolder,extensionstr,tracestr,parameterstr,infermode,algostr)
    if not os.path.exists(tresultfolder):
       os.makedirs(tresultfolder)
    indices = set([-1] + [int(myfile.replace(".hist","")) for myfile in myutil.listfiles(tresultfolder)])
    if complete and len(indices) >= completeupto + 1:
       return
    resultfile = "{0}/{1}.hist".format(tresultfolder,max(indices)+1)
    trunfolder = "{0}/{1}/{2}/{3}/{4}/{5}/{6}".format(runfolder,extensionstr,tracestr,parameterstr,infermode,algostr,max(indices)+1)
    if not os.path.exists(trunfolder):
       os.makedirs(trunfolder)
    tnewtracefolder = "{0}/{1}/{2}/{3}/{4}/{5}/{6}".format(newtracefolder,extensionstr,tracestr,parameterstr,infermode,algostr,max(indices)+1)   
    if not os.path.exists(tnewtracefolder):
       os.makedirs(tnewtracefolder)
    snapshotfile = "{0}/infer.snapshot".format(tnewtracefolder)
    G = InputOutput.readGraphAndParams(path)
    
    if infermode == "History":
       fracpoints = [(timefrac*index)/timecount for index in xrange(1,timecount+1)]
    elif infermode == "Spreader":
       fracpoints = [timefrac]
       for index in xrange(1,timecount):
           interval = (1.0-timefrac)/(timecount-1)
           fracpoints.append(timefrac+(index*interval))
       assert max(fracpoints) <= 1.00000001
    trace = InputOutput.readPlainTrace(tracefile,prob)
    maxtime = Trace.getMaxTraceTime(trace)
    obstimes = sorted(list(set([int(round(frac*maxtime)) for frac in fracpoints])))
    if 0 in obstimes:
       obstimes.remove(0)
    if len(obstimes) == 0:
       return    
    print "observed times" 
    print obstimes
    #if max(obstimes) < 3:
    #   return 
    curstates = [Trace.trace2Snapshot(trace,obstime,smodel,G) for obstime in obstimes]  
    InputOutput.writeSnapshot(curstates,infermode,inter,smodel,obstimes,snapshotfile)
    
    vararr = {"graphfile": path, "dist": prob, "runfolder":trunfolder}
    if inter == None:
       vararr["scoretimes"] = " ".join([str(time) for time in obstimes])
    for item in ["snapshotfile","noise","noisetype","smodel","algo","inter","infermode","evol","printscore","resultfile","tracefile"]:
        exec('vararr["{0}"] = {0}'.format(item)) in locals(),globals()
    for param in algoparams.keys():
        vararr[param] = algoparams[param]
    configpath = "{0}/config_{1}_{2}_{3}_{4}_{5}-{6}.config".format(configfolder,extensionstr,infermode,algostr,parameterstr,tracestr,max(indices)+1)
    genHistoryConfigFile(configpath,vararr)
    code = "python HistoryInfer.py {0}".format(configpath)
    os.system(code)
    return
    if random.random() <= 0.5:
       pool = "pool2"
    else:
       pool = "pool2"     
    pbsfilename = "{0}-{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}-{9}-{10}.pbs".format(realdata,evol,smodel,prob,inter,extensionstr,infermode,algostr,parameterstr,tracestr,max(indices)+1)
    Pbs.submitPbs(code,pbsfolder,pbsfilename,pool)


def assignParams():
    """assign params
    Args:
    Returns:
       paramdict: 
    """
    paramdict = {}
    paramdict["newtracepref"] = "tracesnapshots"
    paramdict["tracepref"] = "traces"
    paramdict["configpref"] = "histconfig"
    paramdict["runpref"] = "run"
    paramdict["resultpref"] = "result"
    paramdict["resultfolderpref"] = "."
    paramdict["runfolderpref"] = "."
    paramdict["tracefolderpref"] = "."
    paramdict["newtracefolderpref"] = "."
    paramdict["graphfolderpref"] = "."
    paramdict["configfolderpref"] = "."
    paramdict["pbsfolder"] = "pbsfolder"
    paramdict["realdata"] = "syn"
    paramdict["evol"] = "static"
    paramdict["infermode"] = "History" #"History" #"Spreader"
    paramdict["smodel"] = "si" #"seir","sir","sis", "samd", "mamd"
    paramdict["prob"] = "dis" #"dis"
    paramdict["inter"] = None #"bound"
    paramdict["startcounts"] = [1,2,3,5,7,10] #[1,2,3,4,5,6,7,8,9,10] #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    paramdict["maxtimefracs"] = [0.15,0.26,0.34,0.4,0.51,0.6,0.76,0.9]
    paramdict["mintimefracs"] = [0.1,0.15,0.2,0.26,0.34,0.4,0.45,0.51,0.55,0.6,0.65,0.7,0.76,0.8,0.85,0.9,0.95] 
    paramdict["timecounts"] = [3,5] #[1,2,3,4,5,6,7] #[1,2,3,4,5,10] #[1,2,3,4,5] #[1,2,3,4,5,10] #[1,2,3,4,5] #[1,2,3,4,5,6,7,10] 
    paramdict["noises"] = [0.0] #[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] 
    paramdict["noisetype"] = "StateChange" #Normal noise
    paramdict["filesamplecount"] = 10 #for syndPlainTrace(tracefile,prob)
    paramdict["complete"] = True
    paramdict["completeupto"] = 2
    paramdict["fraccons"] = {Trace.INFECTED: {"min":0.00001, "max":1.01}}
    paramdict["printscore"] = None #"KenTauCorrelation" #"Graph,rand" #"KenTauCorrelation" #"KenTauCorrelation" #"Hausdorff"
    return paramdict
     

if __name__ == "__main__":
    paramdict = assignParams()
    for var in paramdict.keys():
        if type(paramdict[var]) == type(""):
           exec('{0}="{1}"'.format(var,paramdict[var]))
        else:
           exec('{0}={1}'.format(var,paramdict[var]))               
    graphtraceinput = [graphfolderpref,tracepref,tracefolderpref,newtracepref,newtracefolderpref]
    runresultinput = [runpref,runfolderpref,resultpref,resultfolderpref]
    configinput = [configpref,configfolderpref,pbsfolder]
    sideinput = [smodel,prob,evol,realdata,inter]
    (graphfolder,tracefolder,newtracefolder,configfolder,resultfolder,runfolder) = genMainFolders(graphtraceinput,runresultinput,configinput,sideinput)
    mainparamlist = list(itertools.product(startcounts,noises))
    random.shuffle(mainparamlist)
    for startcount,noise in mainparamlist:
        samplenoisestr = "{0}-{1}-{2}".format(noisetype,noise,0)
        graphinput = [graphfolder,realdata,evol,prob,smodel,filesamplecount,inter]
        traceinput = [tracefolder,fraccons,samplenoisestr,startcount]
        path2info = returnPathInfo(graphinput,traceinput)
        algoblocks = getAlgoBlocks(prob,inter,infermode,smodel)
        if infermode == "Spreader":
           paramlist = list(itertools.product(algoblocks,mintimefracs,timecounts))
        elif infermode == "History":
           paramlist = list(itertools.product(algoblocks,maxtimefracs,timecounts))     
        for path in path2info.keys():
            Gname,tracefile = path2info[path]
            tracestr = "-".join(tracefile.split("/")[-3:])
            for algoblock,timefrac,timecount in paramlist:
                algo,algoparams = algoblock
                if algo in ["NetSleuth","RumorCentrality","KEffectors"]:
                   assert infermode == "Spreader" and prob == "dis"
                   if timecount > 1:
                      continue
                   algofields = [algo,algoparams,runfolder]
                   tracefields = [noise,noisetype,timefrac,timecount,realdata,tracefile,tracestr,newtracefolder,tracefolder]
                   otherfields = [complete,completeupto,path,graphfolder,configfolder,resultfolder]    
                   runOtherAlgos(algofields,tracefields,otherfields) 
                else:
                   algofields = [algo,algoparams,infermode,inter,runfolder]
                   tracefields = [noise,noisetype,timefrac,timecount,prob,smodel,evol,realdata,tracefile,tracestr,newtracefolder,tracefolder]
                   otherfields = [complete,completeupto,path,printscore,graphfolder,configfolder,resultfolder]
                   runMyAlgos(algofields,tracefields,otherfields)
