#  
# This code generates traces for SEIR models
#
#THIS IS NOT SAME AS Cormin TRACE GEN
#
import networkx as nx
import random
import os
import sys
sys.path.append("./lib")
import math
import myutilities as myutil
import cPickle
import gzip
import itertools
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace
from Pbs import Pbs

def genTraceConfigFile(configfile,vararr):
    """generates trace config file
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

def returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount,prob,smodel,graphformat): 
    """returns path info
    Args:
       graphfolder:
       tracefolder:
       ..
    Returns:
       path2info:   
    """      
    path2info = {}
    if realdata == "real" and evol == "static":
       for filename in myutil.listfiles(graphfolder):
           if not filename.endswith(graphformat):
              continue 
           filepath = "{0}/{1}".format(graphfolder,filename)
           if filename.split("-")[0:2] != [prob,smodel]:
              continue
           filetracefolder = "{0}/{1}".format(tracefolder,filename)
           path2info[filepath] = (filename,filetracefolder)       
    elif realdata == "syn" and evol == "static":
       for filename in myutil.listfiles(graphfolder):
           if not filename.endswith(graphformat):
              continue 
           filepath = "{0}/{1}".format(graphfolder,filename)
           if filename.split("-")[0:2] != [prob,smodel]:
              continue
           innum = int(filename.split("_")[-1].replace(".{0}".format(graphformat),""))
           if innum > 1:
              continue 
           filetracefolder = "{0}/{1}".format(tracefolder,filename)
           path2info[filepath] = (filename,filetracefolder)
    elif realdata == "real" and evol == "dynamic":
       for dirname in myutil.listdirectories(graphfolder):
           dirpath = "{0}/{1}".format(graphfolder,dirname)
           if dirname.split("-")[0:2] != [prob,smodel]:
              continue
           dirtracefolder = "{0}/{1}".format(tracefolder,dirname)
           path2info[dirpath] = (dirname,dirtracefolder)
    elif realdata == "syn" and evol == "dynamic":
       for data in myutil.listdirectories(graphfolder):
           dirname = "{0}/{1}".format(graphfolder,data)
           for algo in myutil.listdirectories(dirname):
               algodirname = "{0}/{1}".format(dirname,algo)
               for param in myutil.listdirectories(algodirname):
                   paramdirname = "{0}/{1}".format(algodirname,param)
                   for sampledirname in myutil.listdirectories(paramdirname):
                       samplepath = "{0}/{1}".format(paramdirname,sampledirname)
                       if sampledirname.split("-")[0:2] != [prob,smodel]:
                          continue
                       dirtracefolder = "{0}/{1}".format(tracefolder,sampledirname)
                       path2info[samplepath] = (sampledirname,dirtracefolder)
    return path2info

def genMainFolders(graphpref,graphfolderpref,tracepref,tracefolderpref,configpref,configfolderpref,pbsfolder,realdata,smodel,evol,prob):
    graphfolder = "{0}/{1}_{2}_{3}".format(graphfolderpref,realdata,evol,graphpref)
    tracefolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(tracefolderpref,tracepref,realdata,evol,smodel,"edge",prob)
    configfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(configfolderpref,configpref,realdata,evol,smodel,"edge",prob)
    [os.makedirs(folder) for folder in [pbsfolder,configfolder,tracefolder] if not os.path.exists(folder)]
    return graphfolder,tracefolder,configfolder

def main():    
    """automatically generates traces by calling tracegen.py(assumes trace info is given in graphfile)
    """
    graphpref = "graphs"
    tracepref = "traces"
    configpref = "config"
    graphformat = "edgelist"
    realdata = "syn"
    evol = "static"
    noisetype = "StateChange"
    noises = [0.0]#[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]
    samplerate = 0 
    tracefolderpref = "."
    graphfolderpref = "."
    configfolderpref = "."
    pbsfolder = "pbsfolder"
    smodel = "seir" #"seir","sir"
    prob = "dis" #"dis"
    filesamplecount = 20 #for syn data
    startnodecounts = [1,2,3,4,5,6,7,8,9,10] #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    tracecount = 5
    sismaxtime = 100
   
    if prob == "cont" and samplerate == 0:
       assert noise == 0.0 
    (graphfolder,tracefolder,configfolder) = genMainFolders(graphpref,graphfolderpref,tracepref,tracefolderpref,configpref,configfolderpref,pbsfolder,realdata,smodel,evol,prob)
    path2info = returnPathInfo(graphfolder,tracefolder,realdata,evol,filesamplecount,prob,smodel,graphformat)
    for item,startcount,noise in itertools.product(path2info.items(),startnodecounts,noises):
        path = item[0]
        Gname,filetracefolder = item[1]
        extensionstr = "-".join(path.replace("../","").split("/")[1:])
        noisesamplestr = "{0}-{1}-{2}".format(noisetype,noise,samplerate)
        tracefolder = "{0}/{1}/{2}".format(filetracefolder,noisesamplestr,startcount)
        if os.path.exists(tracefolder):
           existfiles = [filename for filename in myutil.listfiles(tracefolder) if filename.find(".plain") != -1]
           if len(existfiles) >= tracecount:
              continue
        print tracefolder   
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
        #os.system(code)
        #continue
        if random.random() <= 0.5:
           POOL = "pool2"
        else:
           POOL = "pool2"  
        pbsfilename="{0}-{1}-{2}-{3}-{4}-{5}-{6}-{7}-{8}.pbs".format(realdata,evol,smodel,prob,extensionstr,tracecount,samplerate,noise,startcount)
        Pbs.submitPbs(code,pbsfolder,pbsfilename,POOL)
   
if __name__ == "__main__":
    main()
