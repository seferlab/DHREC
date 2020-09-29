import networkx as nx
import unittest
import sys
import os
import itertools
import random
sys.path.append('..')
sys.path.append("../lib")
import myutilities as myutil
from InputOutput import InputOutput
from Trace import Trace
import graphgen
    
class GraphGenTest(unittest.TestCase):
    """graphgen test class
    """
    def __init__(self,testname,G,spreadparam,smodel,prob,filename):
       super(GraphGenTest, self).__init__(testname)
       self.G = G
       self.spreadparam = spreadparam
       self.smodel = smodel
       self.prob = prob
       self.filename = filename

    def readGraphTest(self):
        """ reads graph with attributes
        """
        distpart = "-".join(self.filename.split("/")[-1].split("-")[0:-1])
        tsmodel,tprob,distinfo = Trace.folder2SpreadInfo(distpart)
        self.assertEqual(tsmodel,self.smodel)
        self.assertEqual(tprob,self.prob)
        frozen1 = frozenset(distinfo)
        frozen2 = frozenset(self.spreadparam)
        self.assertEqual(frozen1,frozen2)
        for node1,node2,data in self.G.edges(data=True):
            for field in distinfo.keys():
                if type(distinfo[field]) == dict and distinfo[field].has_key("start"):
                   start = distinfo[field]["start"]
                   end = distinfo[field]["end"]
                   if type(start) == tuple:
                      for index in xrange(len(start)):
                          self.assertLessEqual(start[index],self.G[node1][node2][field][1][index])
                          self.assertGreaterEqual(end[index],self.G[node1][node2][field][1][index])
                   else:
                      self.assertLessEqual(start,self.G[node1][node2][field])
                      self.assertGreaterEqual(end,self.G[node1][node2][field])  
                else:
                   if type(distinfo[field]) == tuple:
                      for index in xrange(len(distinfo[field])): 
                          self.assertLess(abs(distinfo[field][index] - self.G[node1][node2][field][1][index]), 0.00001)
                   else:
                      self.assertLess(abs(distinfo[field] - self.G[node1][node2][field]), 0.00001)
                      
def getManyDist(prob,smodel,eachcount):
    sames = [True,False]
    begin,end = 0.5,8.0
    sbegin,send = 0.0,1.0
    if prob == "cont": 
       dists = ["expo","weibull","rayleigh"]
    elif prob == "dis":
       dists = ["expo","weibull","rayleigh"]    
    random.shuffle(dists)
    spreadparams = []
    if smodel == "si":
       for same,s2i,index in itertools.product(sames,dists,range(0,eachcount)):
           if same:
              if s2i in ["expo","rayleigh"]: 
                 var = random.uniform(begin,end)
                 sprob = random.uniform(sbegin,send)
                 sparam = {Trace.S2I: (s2i, (var,)), Trace.SPROB: sprob}
                 spreadparams.append(sparam)
              elif s2i == "weibull":
                 var,var2 = random.uniform(begin,end),random.uniform(begin,end)
                 sprob = random.uniform(sbegin,send)
                 sparam = {Trace.S2I: (s2i, (var,var2)), Trace.SPROB: sprob}
                 spreadparams.append(sparam) 
           else:   
              if s2i in ["expo","rayleigh"]: 
                 var,var2 = sorted([random.uniform(begin,end) for index in xrange(2)])
                 svar,svar2 = sorted([random.uniform(sbegin,send) for index in xrange(2)])
                 sparam = {Trace.S2I: {"dist": s2i, "start": (var,), "end": (var2,)}, Trace.SPROB: {"start": svar, "end":svar2}}
                 spreadparams.append(sparam)
              elif s2i == "weibull":
                 var,var2 = sorted([random.uniform(begin,end) for index in xrange(2)])
                 var3,var4 = sorted([random.uniform(begin,end) for index in xrange(2)])
                 svar,svar2 = sorted([random.uniform(sbegin,send) for index in xrange(2)])
                 sparam = {Trace.S2I: {"dist": s2i, "start": (var,var2), "end": (var3,var4)}, Trace.SPROB: {"start": svar, "end":svar2}}
                 spreadparams.append(sparam)
    return spreadparams
    

def genMainFolders(graphfolderpref,graphpref,realdata,evol):
    graphfolder = "{0}_{1}_{2}".format(realdata,evol,"graphs")
    [os.makedirs(folder) for folder in [graphfolder] if not os.path.exists(folder)]
    return graphfolder
    

def getTests(G,spreadparam,smodel,prob,graphname):
    """ generates tests for graphgen
    """
    suiteFew = unittest.TestSuite()
    suiteFew.addTest(GraphGenTest("readGraphTest",G,spreadparam,smodel,prob,graphname))
    unittest.TextTestRunner(verbosity = 2).run(suiteFew) 
    return suiteFew

def returnPathInfo(realdata,evol,rawpref,rawfolderpref):
    """ returns path info
    """
    rawgraphfolder = "{0}/{1}_{2}_{3}_graphs".format(rawfolderpref,rawpref,realdata,evol)
    if realdata == "real":
       files = ["{0}/{1}".format(rawgraphfolder,filename) for filename in ["stanford.gml"]]
    elif realdata == "syn":
       files = ["{0}/{1}".format(rawgraphfolder,filename) for filename in myutil.listfiles(rawgraphfolder) if filename.endswith(".gml")]
    return files
    
def main():
    """ main function for testing graphgen
    """
    realdatas = ["syn","real"]
    evols = ["static"]
    probs = ["cont","dis"]
    smodels = ["si","si","seir","sis"]
    graphpref = "graphs"
    graphfolderpref = ".."
    rawfolderpref = ".."
    rawpref = "raw"
    outformat = "pickle"
    eachcount = 100
    
    parts = list(itertools.product(realdatas,evols,probs,smodels))
    random.shuffle(parts)
    for (realdata,evol,prob,smodel) in parts:
        for spreadparam in getManyDist(prob,smodel,eachcount):
            graphfolder = genMainFolders(graphfolderpref,graphpref,realdata,evol)
            paths = returnPathInfo(realdata,evol,rawpref,rawfolderpref)
            for path in paths:
                G = InputOutput.readGraphAndParams(path)
                newgraphfile = graphgen.graphgen(G,path,smodel,spreadparam,prob,graphfolder,outformat)
                newG = InputOutput.readGraphAndParams(newgraphfile,outformat)
                getTests(newG,spreadparam,smodel,prob,newgraphfile)
                os.system("rm -rf {0}".format(newgraphfile))

if __name__ == '__main__':    
    main()
