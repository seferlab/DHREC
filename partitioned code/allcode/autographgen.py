import networkx as nx
import os
import sys
sys.path.append("./lib")
import myutilities as myutil
from Trace import Trace
from InputOutput import InputOutput
from Pbs import Pbs

def retdists(realdata,evol,prob,smodel):
    dists = []
    if prob == "dis":
       if smodel == "si":
          dist = {Trace.S2I: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.2,), "end": (1.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.2,), "end": (1.5,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.1,), "end": (0.2,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.1,), "end": (0.5,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.2,), "end": (0.5,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)

          dist = {Trace.S2I: ("expo", (0.1,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (0.2,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (0.5,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (0.7,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (1.0,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (1.5,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (2.0,)), Trace.SPROB: 1.0}
          dists.append(dist)
          
          #dist = {Trace.S2I: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          #dists.append(dist)
          #dist = {Trace.S2I: ("weibull", (0.5,1.5)), Trace.SPROB: 1.0}
          #dists.append(dist)
          #dist = {Trace.S2I: {"dist": "weibull", "start": (0.5,0.2), "end": (2.0,0.9)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          #dists.append(dist)
          #dist = {Trace.S2I: ("weibull", (0.5,0.5)), Trace.SPROB: 1.0}
          #dists.append(dist)
          
       elif smodel == "sir":
          dist = {Trace.S2I: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.I2R: {"dist": "expo", "start": (0.8,), "end": (1.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.2,), "end": (1.0,)}, Trace.I2R: {"dist": "expo", "start": (0.8,), "end": (1.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.2,), "end": (1.5,)}, Trace.I2R: {"dist": "expo", "start": (0.8,), "end": (1.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
           
          dist = {Trace.S2I: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.I2R: {"dist": "expo", "start": (1.0,), "end": (1.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.2,), "end": (1.0,)}, Trace.I2R: {"dist": "expo", "start": (1.0,), "end": (1.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "expo", "start": (0.2,), "end": (1.5,)}, Trace.I2R: {"dist": "expo", "start": (1.0,), "end": (1.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          
          dist = {Trace.S2I: ("expo", (0.1,)), Trace.I2R: ("expo", (0.1,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (0.2,)), Trace.I2R: ("expo", (0.2,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (0.5,)), Trace.I2R: ("expo", (0.5,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (0.7,)), Trace.I2R: ("expo", (0.5,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (1.0,)), Trace.I2R: ("expo", (0.5,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (1.5,)), Trace.I2R: ("expo", (0.5,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (2.0,)), Trace.I2R: ("expo", (1.0,)), Trace.SPROB: 1.0}
          dists.append(dist)

          dist = {Trace.S2I: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)},  Trace.I2R: ("expo", (0.2,)), Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: ("weibull", (0.5,1.5)), Trace.I2R: ("expo", (0.2,)), Trace.SPROB: 1.0}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)},  Trace.I2R: ("expo", (0.5,)), Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: ("weibull", (0.5,1.5)), Trace.I2R: ("expo", (0.5,)), Trace.SPROB: 1.0}
          dists.append(dist)

          dist = {Trace.S2I: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)},  Trace.I2R: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)},  Trace.I2R: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)}, Trace.SPROB: {"start": 0.8, "end":0.8}}
          dists.append(dist)
          
       elif smodel == "seir":
          dist = {Trace.S2E: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.E2I: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.I2R: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.SPROB: {"start": 1.0, "end":1.0}}
          dists.append(dist)
          
          #dist = {Trace.S2E: ("expo", (0.5,)), Trace.I2R: ("expo", (0.5,)), Trace.E2I: ("expo", (0.5,)), Trace.SPROB: 1.0}
          #dists.append(dist)
          #dist = {Trace.S2E: ("expo", (0.7,)), Trace.I2R: ("expo", (0.5,)), Trace.E2I: ("expo", (0.5,)), Trace.SPROB: 1.0}
          #dists.append(dist)
          #dist = {Trace.S2E: ("expo", (1.0,)), Trace.I2R: ("expo", (0.5,)), Trace.E2I: ("expo", (0.5,)), Trace.SPROB: 1.0}
          #dists.append(dist)
          #dist = {Trace.S2E: ("expo", (1.5,)), Trace.I2R: ("expo", (0.5,)), Trace.E2I: ("expo", (0.5,)), Trace.SPROB: 1.0}
          #dists.append(dist)
          #dist = {Trace.S2E: ("expo", (2.0,)), Trace.I2R: ("expo", (1.0,)), Trace.E2I: ("expo", (1.0,)), Trace.SPROB: 1.0}
          #dists.append(dist)
    elif prob == "cont":
       if smodel == "si":
          dist = {Trace.S2I: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.SPROB: {"start": 0.05, "end":0.3}}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (0.5,)), Trace.SPROB: 0.1}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)}, Trace.SPROB: {"start": 0.05, "end":0.3}}
          dists.append(dist)
          dist = {Trace.S2I: ("weibull", (0.5,1.5)), Trace.SPROB: 0.1}
          dists.append(dist)
       elif smodel == "sir":
          dist = {Trace.S2I: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.I2R: {"dist": "expo", "start": (0.5,), "end": (2.0,)}, Trace.SPROB: {"start": 0.05, "end":0.3}}
          dists.append(dist)
          dist = {Trace.S2I: ("expo", (0.5,)), Trace.I2R: ("expo", (0.5,)), Trace.SPROB: 0.1}
          dists.append(dist)
          dist = {Trace.S2I: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)}, Trace.I2R: {"dist": "weibull", "start": (0.5,1.0), "end": (2.0,2.2)}, Trace.SPROB: {"start": 0.05, "end":0.3}}
          dists.append(dist)
          dist = {Trace.S2I: ("weibull", (0.5,1.5)), Trace.SPROB: 0.1, Trace.I2R: ("weibull", (0.5,1.5))}
          dists.append(dist)
    return dists    

def genGraphGenConfigFile(configfile,vararr):
    """generates graphgen config file
    Args:
       configfile:
       vararr:
    Returns:
    """
    with open(configfile,"w") as file:
        for var in vararr.keys():
            if var == "dists":
               for key in vararr[var].keys():
                   if key in [Trace.S2I, Trace.I2R, Trace.E2I, Trace.I2S, Trace.S2E]:
                      if type(vararr[var][key]) == dict:
                         startstr = " ".join([str(item) for item in vararr[var][key]["start"]])
                         endstr = " ".join([str(item) for item in vararr[var][key]["end"]])
                         file.write("{0} dist: {1}\n".format(key,vararr[var][key]["dist"]))
                         file.write("{0} start: {1}\n".format(key,startstr))
                         file.write("{0} end: {1}\n".format(key,endstr))
                      else:
                         startendstr = " ".join([str(item) for item in vararr[var][key][1]])
                         file.write("{0} dist: {1}\n".format(key,vararr[var][key][0]))
                         file.write("{0} start: {1}\n".format(key,startendstr))
                         file.write("{0} end: {1}\n".format(key,startendstr))
                   elif key == Trace.SPROB:
                      if type(vararr[var][key]) == dict:
                         file.write("{0} start: {1}\n".format(key,vararr[var][key]["start"]))
                         file.write("{0} end: {1}\n".format(key,vararr[var][key]["end"]))
                      else:   
                         file.write("{0} start: {1}\n".format(key,vararr[var][key]))
                         file.write("{0} end: {1}\n".format(key,vararr[var][key]))
            else:
               file.write("{0}: {1}\n".format(var,vararr[var]))


def main():
   pbsfolder = "pbsfolder"
   configfolder = "graphgen_config"
   outformat = "edgelist" 
   realdata = "real" 
   evol = "static"
   prob = "dis" #"dis"
   smodel = "si"
   graphpref = "graphs"
   graphfolderpref = "."
   rawfolderpref = "."
   rawpref = "raw"
   
   if not os.path.exists(pbsfolder):
      os.makedirs(pbsfolder)
   if not os.path.exists(configfolder):
      os.makedirs(configfolder)   
   rawgraphfolder = "{0}/{1}_{2}_{3}_graphs".format(rawfolderpref,rawpref,realdata,evol)
   if realdata == "real":
      files = ["{0}/{1}".format(rawgraphfolder,filename) for filename in myutil.listfiles(rawgraphfolder) if filename.endswith(".gml")]
   elif realdata == "syn":
      files = ["{0}/{1}".format(rawgraphfolder,filename) for filename in myutil.listfiles(rawgraphfolder) if filename.endswith(".gml")]
   if realdata == "syn":
      files = [filename for filename in files if filename.find("_1.gml") != -1 or filename.find("_1.gml") != -1]
   graphfolder = "{0}/{1}_{2}_{3}".format(graphfolderpref,realdata,evol,"graphs")
   alldists = retdists(realdata,evol,prob,smodel)
   
   for gfile in files:
       for dist in alldists:
           vararr = {}
           vararr["graphfilename"] = gfile
           vararr["graphfolder"] = graphfolder
           vararr["smodel"] = smodel
           vararr["outformat"] = outformat
           vararr["dists"] = dist
           vararr["dist"] = prob
           diststr = Trace.getSpreadFolder(smodel,dist,prob)
           gstr = "-".join(gfile.split("/"))
           configpath = "{0}/{1}-{2}-{3}-{4}-{5}.config".format(configfolder,gstr,smodel,outformat,diststr,prob)
           genGraphGenConfigFile(configpath,vararr)
           code = "python graphgen.py {0}".format(configpath)
           #os.system(code)
           #continue
           pbsfilename="{0}-{1}-{2}-{3}-{4}.pbs".format(gstr,smodel,outformat,diststr,prob)
           Pbs.submitPbs(code,pbsfolder,pbsfilename,"pool2")
        
if __name__ == "__main__":
   main()

