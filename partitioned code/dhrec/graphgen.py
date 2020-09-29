import networkx as nx
import os
import sys
sys.path.append("./lib")
import myutilities as myutil
from Trace import Trace
from InputOutput import InputOutput

def assignDefaultGraphGenParameters():
    """assigns default parameter values for graph gen
    Returns:
       returns default value array
    """
    vararr={}
    vararr["graphfilename"] = "graph.gml"
    vararr["graphfolder"] = "graphfolder"
    vararr["smodel"] = "si"
    vararr["outformat"] = "edgelist"
    vararr["dists"] = {}
    vararr["dist"] = "dis"
    return vararr

def readGraphGenConfigFile(filename):
    """reads graphgen configuration file
    Args:
       filename: configuration filename
    Returns:
       returns parameters from configuration file
    """
    vararr = assignDefaultGraphGenParameters()
    with open(filename,"r") as file:
        for line in file:
            varname,value = line.rstrip().split(":")
            value = value.lstrip()
            if len(varname.split(" ")) > 1:
               info, sideinfo = varname.split(" ")
               if info in [Trace.S2I, Trace.I2R, Trace.E2I, Trace.I2S, Trace.S2E]:
                  if not vararr["dists"].has_key(info):
                     vararr["dists"][info] = {}
                  if sideinfo in ["start","end"]:    
                     param = [float(item) for item in value.split(" ")[0:]]
                     vararr["dists"][info][sideinfo] = tuple(param)
                  elif sideinfo == "dist":
                     vararr["dists"][info][sideinfo] = value 
               elif info == Trace.SPROB:
                  if not vararr["dists"].has_key(info):
                     vararr["dists"][info] = {} 
                  vararr["dists"][info][sideinfo] = float(value)
            else:
               vararr[varname] = value
    assert vararr["smodel"] in ["si","sir","seir","sis"] 
    return vararr
   
def graphgen(G,gfile,smodel,dist,prob,graphfolder,outformat):
    """generate graphs with additional features(adds features) in outformat
    Args:
       G: graph 
       gfile: rawfilename
       dist: dist info
       prob: discrete or continuous
       graphfolder: graphfolder
       outformat: graph output format
    Returns: 
       graphfile: new graphfile path
    """
    filename = "{0}-{1}".format(Trace.getSpreadFolder(smodel,dist,prob),gfile.split("/")[-1].replace(".gml",".{0}".format(outformat)))
    graphfile = "{0}/{1}".format(graphfolder,filename)
    if os.path.exists(graphfile):
       return None
    keepkeys = {"si":[Trace.SPROB,Trace.S2I], "sir":[Trace.SPROB,Trace.S2I,Trace.I2R], "seir":[Trace.SPROB,Trace.S2E,Trace.I2R]}
    for key in dist.keys():
        if key not in keepkeys[smodel]:
           del dist[key]        
    InputOutput.WriteGraphAndParams(G,dist,graphfile,smodel,prob,outformat)
    return graphfile
 
def main():
   """adds dist info to graph by the parameters provided in config file
   Args:
      configfile: configuration file
   """
   if len(sys.argv)!=2:
      print "Graphgen Configuration file must be given\n Look at graphgen.config for format details!!"
      exit(1)   
   vararr = readGraphGenConfigFile(sys.argv[1])
   for var in vararr.keys():
       if type(vararr[var])==type(""):
          exec('{0}="{1}"'.format(var,vararr[var]))
       else:
          exec('{0}={1}'.format(var,vararr[var]))             
   assert dist in ["dis","cont"]
   G = InputOutput.readGraphAndParams(graphfilename)
   if not os.path.exists(graphfolder):
      os.makedirs(graphfolder)   
   graphgen(G,graphfilename,smodel,dists,dist,graphfolder,outformat)
        
if __name__ == "__main__":
   main()

