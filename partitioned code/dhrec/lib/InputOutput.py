import networkx as nx
import numpy as np
import scipy as sp
import random
import math
import sys
import os
import myutilities as myutil
import gzip
from copy import deepcopy
import string
import cPickle
import checkUtil
from Trace import Trace


class InputOutput():     
    """ Methods related to Input/Output(graph/trace)
    """
    def __init__(self):
       return
    
    @staticmethod 
    def writeConfigFile(vararr,filename): 
       """write given varar to configuration file
       Args:
          vararr: config file parameters
          filename: configuration filename
       """
       assert vararr["smodel"] in ["si","sir","seir","sis"] 
       with open(filename,"w") as outfile:
          outfile.write("\n".join(["{0}: {1}".format(varname,value) for varname,value in vararr.items()]) + "\n")

    @staticmethod
    def writeGraph2File(retvalues,resultfile,evol):
        """stores retvalues as edgelist to file
        Args:
           retvalues: edge values in dictionary
           resultfilename: graph resultfilename
           evol: static/dynamic graph
        """
        if evol == "static":
           with open(resultfile,"w") as file:
              file.write("".join(["{0} {1} {2}\n".format(node1,node2,retvalues[(node1,node2)]) for node1,node2 in retvalues.keys()]))
        elif evol == "dynamic":
           alltimes = {}
           for node1,node2,time in retvalues.keys():
               alltimes.setdefault(time,set())
               alltimes[time].add((node1,node2))   
           with open(resultfile,"w") as file:
               for time in alltimes.keys():
                   file.write("Time: {0}\n".format(time))
                   file.write("".join(["{0} {1} {2}\n".format(node1,node2,retvalues[(node1,node2)]) for node1,node2 in retvalues[time].keys()]))
               
    @classmethod
    def convertCeferOut(self,cplexoutpath,evol):
        """reads output solution and returns the edges with values
        Args:
           cplexoutpath: CPLEX output file
           evol: static/dynamic graph
        Returns:
           retvalues2: edge variables returned
        """
        retvalues = self.readCplexOut(cplexoutpath,specific = ["x"])[0]
        retvalues2 = {}
        if evol == "dynamic":
           for varname in retvalues.keys():
               node1,node2,time = [int(item) for item in varname.replace("x","").split("?")]
               retvalues2[(node1,node2,time)] = retvalues[varname]
        elif evol == "static":
           for varname in retvalues.keys():
              node1,node2 = [int(item) for item in varname.replace("x","").split("?")]
              retvalues2[(node1,node2)] = retvalues[varname]
        return retvalues2

    @classmethod
    def readCplexOut(self,outfile,specific=[]):
       """reads CPLEX output file and returns only SPECIFIC variable values as dictionary and objective value
        Args:
           outfile: CPLEX output file
           specific: specific variable prefixes such as x
        Returns:
           retvalues: variable-value dictionary
       """
       retvalues = {}
       varflag = False
       objval = None
       with open(outfile,"r") as file:
          for line in file:
              line = line.rstrip()
              if not varflag:
                 if line.startswith("CPLEX>") and line.find("Objective =") != -1 and line.find("Optimal:") != -1:
                    objval = float(line.split("Objective =")[1])
                    continue  
              if not varflag and line.find("CPLEX> Variable Name")!=-1 and line.find("Solution Value")!=-1:
                 varflag=True
                 continue
              if varflag:
                 for varname in specific: 
                     if line.startswith(varname):
                        key,value = line.split()
                        retvalues[key] = float(value)
                        break
       assert objval != None              
       return retvalues,objval
 
    @classmethod
    def addAttribute2Graph(self,G,dists):
        """ appends dist info to G as attributes(if there is interval, assigns one)
        Args:
          G: graph
          dists: distributions
        """
        for key in dists.keys():
            if key in [Trace.I2R,Trace.E2I,Trace.I2S]:
               if type(dists[key]) == dict:
                  for v in G.nodes():
                      G.node[v][key] = (dists[key]["dist"],tuple([random.uniform(dists[key]["start"][index],dists[key]["end"][index]) for index in xrange(len(dists[key]["start"]))]))
               else:
                  for v in G.nodes():
                      G.node[v][key] = dists[key]
            elif key in [Trace.S2I,Trace.S2E]:
               if type(dists[key]) == dict:
                  for u,v in G.edges():
                      G[u][v][key] = (dists[key]["dist"],tuple([random.uniform(dists[key]["start"][index],dists[key]["end"][index]) for index in xrange(len(dists[key]["start"]))]))
               else:
                  for u,v in G.edges():
                      G[u][v][key] = dists[key]
            elif key == Trace.SPROB:
               if type(dists[key]) == dict:
                  for u,v in G.edges():
                      G[u][v][key] = random.uniform(dists[key]["start"],dists[key]["end"]) 
               else:
                  for u,v in G.edges():
                      G[u][v][key] = dists[key]

                      #newdict = {"dist":dists[key]["dist"]}
                      #for index in xrange(len(dists[key]["start"])):
                      #    newdict[index] = random.uniform(dists[key]["start"][index],dists[key]["end"][index]) 
                      #G.node[v][key] = newdict
                      ##G.node[v][key] = [dists[key]["dist"]]+[random.uniform(dists[key]["start"][index],dists[key]["end"][index]) for index in xrange(len(dists[key]["start"]))]
                      #newdict = {"dist":dists[key]["dist"]}
                      #for index in xrange(len(dists[key]["start"])):
                      #    newdict[index] = random.uniform(dists[key]["start"][index],dists[key]["end"][index]) 
                      #G[u][v][key] = newdict
                      ##G[u][v][key] = [dists[key]["dist"]]+ [random.uniform(dists[key]["start"][index],dists[key]["end"][index]) for index in xrange(len(dists[key]["start"]))]
        
    @staticmethod
    def roundCeferOut(edge2val):
        """Reads Graph combined with parameters
        Args:
          edge2val: edges to values
        Returns:
          inferedges: set of inferred edges
        """
        return set([edge for edge in edge2vale.keys() if edge2val[edge] >= 0.000001])
    
    @staticmethod
    def readCeferOut(outfile,evol):
        """Reads CEFER out
        Args:
          outfile: CEFER outfile
          evol: evol
        Returns:
          resdict: result hash
        """
        resdict = {}
        with open(outfile,"r") as infile:
           for line in infile:
               node1,node2,val = line.rstrip().split(" ")
               resdict[(int(node1), int(node2))] = float(val)
        return resdict

    @staticmethod
    def writeCeferOut(outfile,resdict):
        with open(outfile,"w") as outfile:
             outfile.write("\n".join(["{0} {1} {2}".format(node1,node2,val) for (node1,node2),val in resdict.items()]) + "\n")


    @staticmethod
    def readSnapshot(snapshotfile,disttype):
        """reads snapshot from file
           Args:
             snapshotfile:
             disttype: dist/cont
           Returns:
             curstates:
             obstimes:
        """
        curstates,index,smodel,obstimes,obstime,curstate = [], 0, None, [], None,{}
        with open(snapshotfile,"r") as infile:
            for line in infile:
                line = line.rstrip()
                if index == 0:
                   smodel = line
                   index += 1
                   assert smodel in ["si","sir","seir","sis"]
                else:
                   if disttype == "dis":
                      try:
                         val = int(line) 
                         if obstime != None:
                            curstates.append(deepcopy(curstate))
                         obstime = val   
                         obstimes.append(obstime)
                         curstate = {}
                         continue
                      except ValueError:
                         pass
                   elif disttype == "cont":
                      try:
                         val = float(line)  
                         if obstime != None:
                            curstates.append(deepcopy(curstate)) 
                         obstime = val
                         obstimes.append(obstime)
                         curstate = {}
                         continue
                      except ValueError:
                         pass
                   for state in [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]:
                       if line.startswith(state):
                          splitted = line.split("\t")
                          if len(splitted) == 1 :
                             nodeset = set()
                          else:
                             nodeset = set(int(node) for node in splitted[1].split(" ")) 
                          curstate[splitted[0]] = nodeset
            curstates.append(deepcopy(curstate))
            return curstates,obstimes
        
    @staticmethod
    def writeSnapshot(curstates,infermode,inter,smodel,obstimes,snapshotfile):
        """writes snapshot to file
          Args:
            curstates: depends on infermode
            infermode:
            inter:
            smodel:
            obstimes:
            snapshotfile:
        """
        minmap = {None: min(obstimes),"None":min(obstimes),"bound":0}
        with open(snapshotfile,"w") as outfile:
            outfile.write("{0}\n".format(smodel))
            for index in xrange(len(obstimes)):
                outfile.write("{0}\n".format(obstimes[index]-minmap[inter]))
                curstate = curstates[index]
                for state in curstate.keys():
                    nodestr = ""
                    if len(curstate[state]) > 0:
                       nodestr = " ".join([str(node) for node in curstate[state]]) 
                    outfile.write("{0}\t{1}\n".format(state,nodestr))  

    @staticmethod                  
    def history2Times(sol,seenstates,smodel):
        """converts inferred to hist to state transition times
        Args:
           sol:
           seenstates:
           smodel:
        Returns:
           newsol: dict showing time and nodes transitioned at that times
        """
        newsol = {time:{state:set() for state in  [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]} for time in sol.keys()}
        allnodes = set(node for time in sol.keys() for state in sol[time].keys() for node in sol[time][state])
        seentimes = {node: {} for node in allnodes}
        for time in sorted(sol.keys()):
            for state in [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]:
                for node in sol[time][state]:
                    if not seentimes[node].has_key(state):
                       seentimes[node][state] = time 
        for node in seentimes.keys():
            for state in seentimes[node].keys(): 
                newsol[seentimes[node][state]][state].add(node)
        return newsol
                   
    @staticmethod                  
    def history2TimesNOTUSED(sol,seenstate,smodel):
        """converts inferred to hist to state transition times
        Args:
           sol:
           seenstate:
           smodel:
        Returns:
           newsol: dict showing time and nodes transitioned at that times
        """
        sol[0] = deepcopy(seenstate)
        allnodes = set(node for time in sol.keys() for state in sol[time].keys() for node in sol[time][state])
        soltimes = sorted(sol.keys())
        transtimes = {Trace.EXPOSED:{},Trace.INFECTED:{},Trace.RECOVERED:{}}
        for time in soltimes:
            for state in sol[time].keys():
                if state == Trace.SUSCEPTIBLE:
                   continue
                for node in sol[time][state]:
                    transtimes[state].setdefault(node,set()).add(time)
        print "in his2times sol"
        for time in sorted(sol.keys()):
            for state in sol[time].keys(): 
                print time,state,sol[time][state]           
        print "hist2times"            
        for state in transtimes.keys():
            for node in transtimes[state].keys():
                print state,node,transtimes[state][node]
                transtimes[state][node] = min(transtimes[state][node])
        newsol = {time: {state: set() for state in [Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]} for time in soltimes} 
        for state in transtimes.keys():
            for node in transtimes[state].keys():
                newsol[transtimes[state][node]][state].add(node) 
                #newsol[transtimes[state][node]].setdefault(state,set()).add(node) 
        print "newsol"
        for time in newsol.keys():
            for state in newsol[time].keys():
                print time,state,newsol[time][state]          
        return newsol
                
    @staticmethod
    def writeHistoryResult(history,infermode,resultfile,wformat="hist"):
        """writes discrete greedy solution to file 
          Args:
            history: depends on infermode
            infermode:
            resultfile:
            wformat:
        """
        if resultfile.endswith(".pkl"):
           useformat = "pkl"
        elif resultfile.endswith(".hist"):
           useformat = "hist"
        else:
           useformat = "{0}".format(wformat) 
        if useformat == "pkl":
           if infermode == "Spreader":  
              with gzip.open(resultfile,"wb") as outfile:
                 cPickle.dump(history[1])
                 cPickle.dump(history[0])
           elif infermode == "History":
              with gzip.open(resultfile,"wb") as outfile:
                 cPickle.dump(history) 
        elif useformat == "hist":
           if infermode == "Spreader":
              with open(resultfile,"w") as outfile:
                 outfile.write("{0}\n".format(history[1]))  
                 for state in history[0].keys():
                     outfile.write("{0}\t{1}\n".format(state," ".join([str(node) for node in history[0][state]])))
           elif infermode == "History":
              with open(resultfile,"w") as outfile:
                 for time in sorted(history.keys()):
                     for state in history[time].keys():
                         infostr = ""
                         if len(history[time][state]) != 0:
                            infostr = " ".join([str(node) for node in history[time][state]])
                         outfile.write("{0}\t{1}\t{2}\n".format(time,state,infostr))
                        
    @staticmethod
    def readHistoryResult(resultfile,infermode,rformat="hist"):
        """reads discrete greedy output from file 
          Args:
            resultfile: resultfile to be read
            infermode:
            rformat: read format
          Returns:
            history: depends on infermode
        """
        if resultfile.endswith(".pkl"):
           useformat = "pkl"
        elif resultfile.endswith(".hist"):
           useformat = "hist"
        else:
           useformat = "{0}".format(rformat)
        if useformat == "pkl":   
           if infermode == "Spreader": 
              with open(resultfile,"rb") as infile:
                 time = cPickle.load(infile)
                 instate = cPickle.load(infile)
              return (time,instate)   
           elif infermode == "History": 
              with open(resultfile,"rb") as infile:
                 history = cPickle.load(infile)
              return history
        elif useformat == "hist":
           if infermode == "Spreader":
              index,time,instate = 0,None,{} 
              with open(resultfile,"r") as infile:
                 for line in infile:
                     line = line.rstrip() 
                     if index == 0:
                        time = int(line)
                        index += 1
                     else:
                        splitted = line.split("\t")
                        if len(splitted) != 1:
                           instate[splitted[0]] = set(int(item) for item in splitted[1].split(" "))
                        else:
                           instate[splitted[0]] = set() 
              return (time,instate)
           elif infermode == "History":
              history = {} 
              with open(resultfile,"r") as infile:
                 for line in infile:
                     line = line.rstrip()
                     splitted = line.split("\t")
                     time,state = int(splitted[0]),splitted[1]
                     if len(splitted) == 3:
                        nodeset = set(int(item) for item in splitted[2].split(" "))
                     else:
                        nodeset = set()
                     history.setdefault(time,{})
                     history[time][state] = nodeset
              return history 
                 
    @staticmethod
    def readGraphAndParams(graphpath,outformat="pickle"):
       """Reads Graph combined with parameters
       Args:
          graphpath: graph file
       Returns:
          G: directed graph with attributes
       """
       if graphpath.endswith("pickle"):
          outformat = "pickle"
       elif graphpath.endswith(".gml"):
          outformat = "gml" 
       elif graphpath.endswith(".edgelist"):
          outformat = "edgelist"    
       if outformat == "pickle": 
          return nx.DiGraph(nx.read_gpickle(graphpath))
       elif outformat == "gml":
          return nx.DiGraph(nx.read_gml(graphpath,relabel=True))
       elif outformat == "edgelist":
          return nx.DiGraph(InputOutput.readGraphEdgeListWithParams(graphpath)) 

    @classmethod   
    def readGraphEdgeListWithParams(self,inpath):
        """read given graph with attr. as edgelist format from inpath
        Args:
           inpath: out graphpath 
        Returns:
           G: graph
        """
        G = nx.DiGraph()
        attr = False
        with open(inpath,"r") as infile:
            for line in infile:
                if len(line.rstrip().split("\t")) == 3:
                   attr = True
                   break
        if attr:        
           with open(inpath,"r") as infile:
              for line in infile:
                  splitted = line.rstrip().split("\t")
                  if len(splitted) == 1: #node
                     G.add_node(int(splitted[0]))
                  elif len(splitted) == 2: #node
                     node,diststr = line.split("\t")
                     newstr = str("")+"-"+str("")+"-"+"_".join(diststr.split(" ")) + "_s2i_expo_0.5"
                     tsmodel,tprob,distinfo = Trace.folder2SpreadInfo(newstr)
                     del distinfo["s2i"]
                     node = int(node)
                     G.add_node(node)
                     G.node[node] = deepcopy(distinfo)
                  elif len(splitted) == 3: #edge   
                     node1,node2,diststr = line.split("\t")
                     newstr = str("")+"-"+str("")+"-"+"_".join(diststr.split(" "))
                     tsmodel,tprob,distinfo = Trace.folder2SpreadInfo(newstr)
                     node1,node2 = int(node1),int(node2)
                     G.add_edge(node1,node2)
                     G[node1][node2] = deepcopy(distinfo)
                  else:
                     print "Edgelist parsing error!!"
                     exit(1)    
        else:
           with open(inpath,"r") as infile:
              for line in infile:
                  splitted = line.rstrip().split("\t")
                  if len(splitted) == 1: #node without attribution
                     G.add_node(int(splitted[0]))
                  elif len(splitted) == 2: #edge without attribution
                     G.add_edge(int(splitted[0]),int(splitted[1]))               
        return G

    @classmethod
    def writeGraphEdgeListWithParams(self,G,outpath,smodel,prob,ongraph):
        """Writes given graph with attributes as edgelist format(if no params write simple graph with also nodes)
        Args:
           G:
           outpath: out graphpath 
           smodel: spreading model
           prob:
           ongraph: attribute exists or not
        """
        with open(outpath,"w") as outfile:
            if ongraph:
               for node in G.nodes():
                   for label in ["id","label"]:
                       if G.node[node].has_key(label):
                          del G.node[node][label]
                   splitted = Trace.getSpreadFolder(smodel,G.node[node],prob).split("_")
                   outstr = " ".join(splitted[0].split("-")[2:3]+splitted[1:])
                   outfile.write("{0}\t{1}\n".format(node,outstr)) 
               for node1,node2 in G.edges():
                   splitted = Trace.getSpreadFolder(smodel,G[node1][node2],prob).split("_")
                   outstr = " ".join(splitted[0].split("-")[2:3]+splitted[1:])
                   outfile.write("{0}\t{1}\t{2}\n".format(node1,node2,outstr))
            else:
               outfile.write("\n".join([str(node) for node in G.nodes()])+"\n")
               outfile.write("\n".join(["{0}\t{1}".format(node1,node2) for node1,node2 in G.edges()])+"\n")
                               
    @classmethod  
    def WriteGraphAndParams(self,G,dists,outfile,smodel,dist,outformat="pickle"):
       """Writes Graph combined with parameters as weighted graph to outfile
       Args:
         G: graph
         dists: model distributions
         outfile: graph filename
         smodel: spreading model
         dist:
         outformat:
       """
       self.addAttribute2Graph(G,dists)
       if outformat == "pickle": 
          nx.write_gpickle(G,outfile)
       elif outformat == "gml":
          nx.write_gml(G,outfile)
       elif outformat == "edgelist":
          self.writeGraphEdgeListWithParams(G,outfile,smodel,dist,True)

    @staticmethod
    def writeTrace(trace,tracefilename,outformat="plain"):
        """writes traces to outfilename
        Args:
           trace: trace
           outfilename: trace outfilename
           smodel:
        """
        methodname = "write{0}Trace".format(outformat.capitalize())
        getattr(InputOutput,methodname)(trace,outfilename)
    
    @staticmethod
    def writePklTrace(trace,outfilename):
       """writes traces to outfilename
       Args:
          trace: trace
          outfilename: trace outfilename
          smodel:
       """
       with gzip.open(outfilename,"wb") as outfile:
           cPickle.dump(trace,outfile)

    @staticmethod
    def writePlainTrace(trace,outfilename):
       """Writes traces to readable text file as plain
       Args:
          trace: trace
          outfilename:
          smodel:
       """
       with open(outfilename,"w") as file:
          if Trace.IsTraceNoisy([trace]):
             sortedtimes = sorted(trace[trace.keys()[0]].keys())
             tracestr = "\n".join(["{0} {1} {2} {3}".format(node,time,state,trace[node][time][state]) for node in trace.keys() for time in sortedtimes for state in trace[node][time].keys()])
          else:
             unsortednodes = []
             for node in trace.keys():
                 if trace[node].has_key(Trace.INFECTED):
                    unsortednodes.append((node,trace[node][Trace.INFECTED]))
                 elif trace[node].has_key(Trace.EXPOSED):
                    unsortednodes.append((node,trace[node][Trace.EXPOSED]))
             sortednodes = [node for node,time in sorted(unsortednodes, key = lambda element : element[1])]
             assert len(set(sortednodes).intersection(set(trace.keys()))) == len(sortednodes)
             tracestr = "\n".join(["{0} {1} {2}".format(node,state,trace[node][state]) for node in sortednodes for state in trace[node].keys()])  
          file.write("{0}\n".format(tracestr))
    
    @staticmethod
    def readPlainTrace(tracefile,traceinfo="cont"):
        """reads given tracefile
        Args:
          tracefile: trace file
          traceinfo: dis/cont
        """
        if traceinfo == "dis":
           method = int
        elif traceinfo == "cont":
           method = float 
        trace = {}
        with open(tracefile,"r") as file:
          for line in file:
              if len(line.rstrip().split("\n")[0].split(" ")) == 3:
                 noisy = False 
              elif len(line.rstrip().split("\n")[0].split(" ")) == 4:
                 noisy = True
        with open(tracefile,"r") as file:
          for line in file:
              for item in line.rstrip().split("\n"):
                  splitted = item.split(" ")
                  if noisy:
                     trace.setdefault(int(splitted[0]),{})
                     trace[int(splitted[0])].setdefault(int(splitted[1]),{})
                     trace[int(splitted[0])][int(splitted[1])][splitted[2]] = method(splitted[3])
                  else:
                     trace.setdefault(int(splitted[0]),{})
                     trace[int(splitted[0])][splitted[1]] = method(splitted[2])
        return trace

    @staticmethod
    def readPklTrace(tracefile):
        """ reads given tracefile
        """
        with gzip.open(tracefile,"rb") as file:
            return cPickle.load(file)

    @staticmethod
    def readTraces(tracefolder,tracecount,smodel,mininfected=-1,informat="plain",rettype="list"):
        """reads tracecount many traces
        Args:
          tracefolder: tracefolder
          tracecount: number of traces wanted
          smodel: spreading model
          mininfected: minimum number of infected nodes in each trace
          informat: trace format
          rettype = return type
        Returns:
          traces: traces as dict
        """
        assert rettype in ["list","dict"] 
        if smodel in ["si","sir","sis"]:
           statekey = Trace.INFECTED
        elif smodel == "seir":
           statekey = Trace.EXPOSED
        allnodes = set()
        if rettype == "list":
           traces = []
        elif rettype == "dict":
           traces = {} 
        filenames = myutil.listfiles(tracefolder)
        random.shuffle(filenames)  
        for filename in filenames:
            tracefilepath = "{0}/{1}".format(tracefolder,filename)
            method = getattr(InputOutput,"read{0}Trace".format(informat.capitalize()))
            trace = method(tracefilepath)
            if len([node for node in trace.keys() if trace[node].has_key(statekey)]) >= mininfected:
               if rettype == "list": 
                  traces.append(trace)
               elif rettype == "dict":
                  traces[filename] = trace 
               allnodes |= set(trace.keys())
            if len(traces) == tracecount:
               return (traces, allnodes) 
        print "Error!! NOT ENOUGH TRACES under {0}".format(tracefolder)
        exit(1)
