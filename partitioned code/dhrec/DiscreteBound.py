import networkx as nx
import numpy as np
import scipy as sp
import scipy.optimize
import random
import math
import sys
sys.path.append("./lib")
import os
import operator
from copy import deepcopy
import myutilities as myutil
from Trace import Trace
from Constraint import Constraint
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from HistoryEnsemble import DiscreteNoneEnsemble as DisEns
from DiscreteDistribution import DiscreteDistribution as Dist
import HistoryUtil
import checkUtil
import boundInferUtil
from Qsapmin import QsapMinTest

class BoundTest():
     #Bound Test Class
     def __init__():
         return

     @staticmethod
     def checkDiffusionValidity(curassign,labelnodeG,mode):
         """checks diffusion validity
         Args:
            curassign:
            labelnodeG:
            mode: label/var
         Returns:
            bool: true/false
         """
         assert mode in ["var","label"]
         if mode == "var":
            nodelabelpairs = [(node,curassign[node].replace("x{0}?".format(node),"")) for node in curassign.keys()]
         elif mode == "label":
            nodelabelpairs = [(node,curassign[node]) for node in curassign.keys()]    
         curG = nx.DiGraph(labelnodeG.subgraph(nodelabelpairs))
         for node,label in nodelabelpairs:
             if len(labelnodeG.predecessors((node,label))) != 0 and len(curG.predecessors((node,label))) == 0:
                return False
         return True  
       
     @staticmethod
     def checkSolution(linear,quad,sol,objval,labelnodeG,solmode):
         """checks solution
         Args:
            linear:
            quad:
            sol:
            objval:
            labelnodeG:
            solmode: var or label
         Returns:
         """
         linsum = sum([-1.0*linear[var] for var in sol.values() if linear.has_key(var)])
         assert linsum <= 0.0 and objval <= 0.0 and linsum >= objval
         assert BoundTest.checkDiffusionValidity(sol,labelnodeG,solmode)
         return True

     @staticmethod
     def checkCoef(linear,quad):
         """checks coefficients
         Args:
            linear:
            quad:
         Returns:
            bool: true/false 
         """
         for var1,var2 in quad.keys():
             assert not quad.has_key((var2,var1))
         for val in quad.values():
             assert val >= 0.0
         for val in linear.values():
             assert val >= 0.0 
         return True
     
class BoundObjFunc():
     #Objective Function class for bounded case
     def __init__(self,G,smodel,curstates,algo):
         """constructor
         Args:
            G:
            smodel:
            curstates: dictionary keyed by obstimes
            algo:
         Returns:
         """
         assert algo in ["search","continuous"]
         self.algo = algo
         self.G = nx.DiGraph(G)
         self.smodel = smodel
         if self.smodel in ["si","sis"]:
            self.states = [Trace.SUSCEPTIBLE,Trace.INFECTED]
         elif self.smodel == "sir":
            self.states = [Trace.SUSCEPTIBLE,Trace.INFECTED,Trace.RECOVERED]
         elif self.smodel == "seir":
            self.states = [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]
         if smodel in ["si","sir","sis"]:
            self.statekey = Trace.INFECTED
            self.transdist = Trace.S2I
         elif smodel == "seir":
            self.statekey = Trace.EXPOSED 
            self.transdist = Trace.S2E    
         self.curstates = deepcopy(curstates)
         self.LOGZERO = -10000000000.0 
         self.ZERO = 0.00000001
         self.EPSILON = 100.0
         self.maxtime = max(curstates.keys())
         node2slots = boundInferUtil.assignStateTimeSlots(curstates,smodel,G)
         self.node2vars = boundInferUtil.getNode2Vars(node2slots,smodel)
         self.labelnodeG = boundInferUtil.makeLabelNodeGraph(self.G,self.node2vars,smodel,"normal","label")
         self.estimateCoefs()
         self.presolution = {} #solution obtained by coloops
         self.findColoops()
         self.var2index,self.index2var = {}, {}
         index = 0
         for node in self.node2vars.keys():
             for var in self.node2vars[node]:
                 self.var2index[var] = index
                 self.index2var[index] = var
                 index += 1
         self.labelcount = index
         self.SAMPLECOUNT = 2 * self.labelcount
         self.MAXNOMINATOR = sum(self.linear.values()) + sum([-1.0*math.log(val) for vals in self.quad.values() for val in vals])
         if globals()["isTest"]:
            BoundTest.checkCoef(self.linear,self.quad) 
            [QsapMinTest.checkVariableValidity(self.curstates,varname,self.smodel) for node in self.node2vars.keys() for varname in self.node2vars[node]]
            assert self.labelcount == sum([len(self.node2vars[node]) for node in self.node2vars.keys()]) and self.labelcount == len(self.index2var.keys())
            assert sum([len(self.node2vars[node]) for node in self.node2vars.keys()]) + len(self.presolution.keys()) == self.labelnodeG.number_of_nodes()
            assert len(set(self.presolution.keys()).intersection(set(self.node2vars.keys()))) == 0
            
     def findColoops(self):
         """eliminates the items that are seen in every basis
         Args:
         Returns:
         """ 
         coloops = [node for node in self.node2vars.keys() if len(self.node2vars[node]) == 1]
         for node in coloops:
             var = list(self.node2vars[node])[0]
             label = var.replace("x{0}?".format(node),"")
             self.presolution[node] = label
             del self.node2vars[node]
             
     def estimateCoefs(self):
         """estimates and stores linear and quads
         Args:
         Returns:
         """
         snodes = set(self.G.nodes()).difference(set(self.node2vars.keys()))
         self.linear = {var: 0.0 for node in self.node2vars.keys() for var in self.node2vars[node]}
         self.quad = {}
         for node in self.node2vars.keys():
             for var1 in self.node2vars[node]:
                 times1 = HistoryUtil.var2Times(var1)
                 for neighnode in set(self.G.predecessors(node)).intersection(set(self.node2vars.keys())):
                     for var2 in self.node2vars[neighnode]:
                         times2 = HistoryUtil.var2Times(var2)
                         if HistoryUtil.checkAffect(times2,times1,self.smodel,Trace.INFECTED,self.statekey,"normal"):
                            f1 = Dist.genPartDist(self.G[neighnode][node][self.transdist][0],self.G[neighnode][node][self.transdist][1],"reverseCdf",self.G[neighnode][node][Trace.SPROB])             
                            fdist = Dist.genPdf(f1)
                            fdist[0] = 1.0
                            diftime = times1[self.statekey] - times2[Trace.INFECTED]
                            assert diftime > 0 
                            first,second = self.ZERO, self.ZERO
                            if fdist.has_key(diftime-1) and fdist[diftime-1] > self.ZERO:
                               first = fdist[diftime-1]
                            if fdist.has_key(diftime) and fdist[diftime] > self.ZERO:
                               second = fdist[diftime]
                            self.quad[(var1,var2)] = (first,second)
                 if self.smodel in ["sir","seir"]:
                    coef = boundInferUtil.estimateEndoCoefs(node,times1,self.G,"i2r",self.ZERO,self.maxtime)
                    if coef != 0.0:
                       self.linear[var1] += coef 
                 if self.smodel == "seir":
                    coef = boundInferUtil.estimateEndoCoefs(node,times1,self.G,"e2i",self.ZERO,self.maxtime)
                    if coef != 0.0:
                       self.linear[var1] += coef 
         locallinear = boundInferUtil.estimateNoneTransCoef(snodes,self.G,self.node2vars,self.smodel,self.ZERO,self.maxtime)
         for key in locallinear.keys():
             self.linear[key] += locallinear[key]
         for key in self.linear.keys():
             assert self.linear[key] >= 0.0
             if self.linear[key] == 0.0 or self.linear[key] == -0.0:
                del self.linear[key] 

     def fracObjValNotUsed(self,fracsol):
         """estimates fractional objective value
         Args:
            fracsol: frac sol
         Returns:
            fracval: fractional objective value
         """
         fracval = 0.0 
         for node in fracsol.keys():
             for var1 in fracsol[node]:
                 firstpart,secondpart = [],[]
                 times1 = HistoryUtil.var2Times(var1)
                 if times1.has_key(Trace.INFECTED) and times1[Trace.INFECTED] == 0:
                    continue
                 if self.smodel in ["si","sir"]:
                    assert times1.has_key(Trace.INFECTED)
                 elif self.smodel == "seir":
                    assert times1.has_key(Trace.EXPOSED) or times1.has_key(Trace.INFECTED)    
                 for neighnode in set(self.G.predecessors(node)).intersection(set(fracsol.keys())):
                     for var2 in fracsol[neighnode]:
                         var2 = fracsol[neighnode]
                         if not self.quad.has_key((var1,var2)):
                            continue 
                         firstpart.append(self.quad[(var1,var2)][0]**fracsol[neighnode][var2])
                         secondpart.append(self.quad[(var1,var2)][1]**fracsol[neighnode][var2])
                 if len(firstpart) == 0:
                    return self.LOGZERO 
                 logval = math.log(reduce(operator.mul,firstpart,1.0) - reduce(operator.mul,secondpart,1.0))
                 assert logval <= 0
                 fracval += (fracsol[node][var1]*logval)
         fracval -= sum([0.0]+[self.linear[var] for node in fracsol.keys() for var in fracsol[node] if self.linear.has_key(var)])
         return fracval
     
     def getBinaryObjVal(self,samplelabel):
         """estimates binary obj func not frac
         Args:
            samplelabel: is of type var(not label)
         Returns:
            binaryval:
         """
         binaryval = 0.0 
         for node in samplelabel.keys():
             var1 = samplelabel[node]
             if var1 == None:
                continue 
             firstpart,secondpart = [],[]
             times1 = HistoryUtil.var2Times(var1)
             if times1.has_key(Trace.INFECTED) and times1[Trace.INFECTED] == 0:
                continue
             if self.smodel in ["si","sir"]:
                assert times1.has_key(Trace.INFECTED)
             elif self.smodel == "seir":
                assert times1.has_key(Trace.EXPOSED) or times1.has_key(Trace.INFECTED)    
             for neighnode in set(self.G.predecessors(node)).intersection(set(samplelabel.keys())):
                 var2 = samplelabel[neighnode]
                 if not self.quad.has_key((var1,var2)):
                    continue 
                 firstpart.append(self.quad[(var1,var2)][0])
                 secondpart.append(self.quad[(var1,var2)][1])
             if len(firstpart) == 0:
                binaryval += self.LOGZERO    
             logval = math.log(reduce(operator.mul,firstpart,1.0) - reduce(operator.mul,secondpart,1.0))
             assert logval <= 0
             binaryval += logval
         binaryval -= sum([0.0]+[self.linear[var] for var in samplelabel.values() if self.linear.has_key(var)])
         return binaryval

     def getRandomLabel(self,var2val):
         """returns randomized var type label
         Args:
            label2val:
         Returns:
            foundvar:
         """
         sumval = float(sum(var2val.values()))
         curvars = var2val.keys()
         random.shuffle(curvars)
         p = random.random() 
         curval = 0.0
         foundvar = None
         for curvar in curvars:
             curval += var2val[curvar]/sumval
             if curval >= p:
                foundvar = curvar
                break  
         if foundvar == None:
            foundvar = curvars[-1]
         return foundvar
          
     def sampleObjVal(self,fracsol):
         """approximates objective function via sampling samplecount times
         Args:
            fracsol:
         Returns:
            estval:
         """   
         estval = 0.0
         for index in xrange(self.SAMPLECOUNT):
             samplelabel = {node : self.getRandomLabel(fracsol[node]) for node in fracsol.keys()}
             estval += self.getBinaryObjVal(samplelabel)  
         estval /= float(self.SAMPLECOUNT)
         assert estval <= 0.0 
         return estval

     def retBasePackingNumber(self):
         """returns matroid based packing number
         Args:
         Returns:
            k:
         """
         k = min([len(self.node2vars[node]) for node in self.node2vars.keys()])
         assert k > 1
         return k

     def isInsideBt(self,cursol,t):
         """checks whether current solution obeys bt constraints
         Args:
            cursol:
            t:
         Returns:
            bool: true/false
         """
         for node in cursol.keys():
             for val in cursol[node].values():
                 if val > t:
                    return False 
         return True

     def checkBaseConstraint(self,initsol,right,mode):
         """checks base constraint of partition matroid(equality =1 )
         Args:
            initsol:
            rho:
            mode: exact or inside
         Returns:
            bool: true/false
         """
         assert mode in ["exact","inside"]
         for node in initsol.keys():
             cursum = sum([initsol[node][label] for label in initsol[node].keys()])
             if mode == "exact":
                if cursum != right:
                   return False  
             elif mode == "inside":
                if cursum > right:
                   return False     
         return True
     
     def findInitialSolution(self,basenum):
         """finds initial solution
         Args:
            basenum:
         Returns:
            initsol: initsol format is dict of dict
         """
         bases = self.findIndependentBases(basenum)
         t = (basenum-1)/float(basenum)
         initsol = {node: {var:0.0 for var in self.node2vars[node]} for node in self.node2vars.keys()}
         for index in bases.keys():
             for var in bases[index]:
                 node = int(var.split("?")[0].replace("x",""))
                 initsol[node][var] = 1 
         if globals()["isTest"]:
            lens = [len(self.node2vars[node]) for node in self.node2vars.keys()]
            minbase,maxbase = min(lens), max(lens)
            seenvals = set(initsol[node][var] for node in initsol for var in initsol[node].keys())
            if minbase != maxbase:
               assert len(seenvals) == 2 and 0.0 in seenvals
            else:
               assert len(seenvals) == 1
            zeronodes = set(node for node in self.node2vars.keys() for var in self.node2vars[node] if initsol[node][var] == 0.0)
            bignodes = set(node for node in self.node2vars.keys() if len(self.node2vars[node]) > min(lens))
            assert len(bignodes.symmetric_difference(zeronodes)) == 0 
            assert self.isInsideBt(initsol,basenum-1) and self.checkBaseConstraint(initsol,basenum,"inside")
         return initsol
         
     def findIndependentBases(self,basenum):
         """finds independent bases
         Args:
            basenum:
         Returns:
            bases:
         """
         bases = {index: [] for index in xrange(basenum)}
         for node in self.node2vars.keys():
             curnodes = list(self.node2vars[node]) 
             random.shuffle(curnodes)
             for index in xrange(basenum):
                 bases[index].append(curnodes[index])
         if globals()["isTest"]:   
            assert basenum <= min([len(self.node2vars[node]) for node in self.node2vars.keys()])      
            for index1 in bases.keys():
                assert len(bases[index]) == len(self.node2vars.keys())
                for index2 in bases.keys():     
                    if index1 == index2:
                       continue
                    assert len(set(bases[index1]).intersection(set(bases[index2]))) == 0 
         return bases
     
     def isDirectionExists(self,cursol,curval,basenum):
         """checks whether improving direction exists
         Args:
            cursol:
            curval: solution of curval
            basenum:
         Returns:
            bool: is new solution exist or not
            newsol:
            newsolval:
         """
         for index1 in xrange(self.labelcount):
             var1 = self.index2var[index1]
             node1 = int(var1.split("?")[0].replace("x",""))
             for index2 in xrange(self.labelcount):
                 var2 = self.index2var[index2]
                 node2 = int(var2.split("?")[0].replace("x",""))
                 if index1 == index2:
                    continue 
                 if cursol[node1][var1] == basenum-1 or cursol[node2][var2] == 0:
                    continue
                 cursol[node1][var1] += 1
                 cursol[node2][var2] -= 1
                 cursolval = self.sampleObjVal(cursol)
                 OPT = self.MAXNOMINATOR / (float(self.labelcount)**2)
                 if (cursolval - curval) >= (1.0/(basenum*self.SAMPLECOUNT)) * OPT: 
                    return True,cursol,cursolval
                 cursol[node1][var1] -= 1
                 cursol[node2][var2] += 1
         return False,cursol,curval

     def contOptimizeObj(self):
         """continuous optimization from Jan Vondrak paper
         Args:
         Returns:
            sol:
         """
         basenum = self.retBasePackingNumber() 
         cursol = self.findInitialSolution(basenum)
         curobjval = self.sampleObjVal(cursol)
         while True:
             existflag,nextsol,nextobjval = self.isDirectionExists(cursol,curobjval,basenum)
             if existflag:
                cursol = deepcopy(nextsol)
                curobjval = nextobjval 
             else:
                break
         if globals()["isTrue"]:    
            assert self.checkBaseConstraint(cursol,basenum,"inside")
         print cursol
         while True:
            solnode2label = self.pipageRound(cursol) 
            for node in self.presolution.keys():
                solnode2label[node] = self.presolution[node]
            if BoundTest.checkDiffusionValidity(solnode2label,self.labelnodeG,"label"):
               break
         sol = {time: {state: set() for state in [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]} for time in xrange(len(self.maxtime+1))}   
         for node in solnode2label.keys():
             varname = "x{0}?".format(solnode2label[node])
             times = HistoryUtil.var2Times(varname)
             for state in times.keys():
                 sol[times[state]][state].add(node) 
         return sol

     def optimizeObj(self):
         """optimizes the objective F
         Args:
         Returns:
            node2label:
         """
         if self.algo == "search":
            return self.searchOptimizeObj()
         elif self.algo in ["continuous","direct"]:
            return self.contOptimizeObj()

     def hitConstraint(self,y,i,j):
         """find hitting constraint
         Args:
            y:
            i:
            j:
         Returns:
            ynew:
            A:
         """ 
         ynew,A = None,None
         for node in self.node2vars.keys():
             if i == 0 and j == 1:
                break 
                              
     def roundPipage(self,fracsol,basenum):
         """does pipage rounding
         Args:
            fracsol:
         Returns:
            node2label:
         """  
         #v,r,q = basenum/float(basenum-1),basenum-1,basenum
         #t = 1/v
         #rho = 1.0/q 
         sol = {}
         node2label = {}
         while len(sol.keys()) != len(self.node2vars.keys()):
            T = list(set(self.node2index.keys()).difference(sol.values()))  
            while len(T) != 0:
                random.shuffle(T)
                i,j = T[0:2]
                yplus,Aplus = self.hitConstraint(y,i,j)
                yminus,Aminus = self.hitConstraint(y,j,i)
                norm1 = np.linalg.norm([yplus[index] - y[index] for index in xrange(len(yplus))])
                norm2 = np.linalg.norm([yplus[index] - yminus[index] for index in xrange(len(yplus))])
                if random.random() < norm1/norm2:
                   y = yminus
                   T = list(set(T).intersection(Aminus))
                else:
                   y = yplus 
                   T = list(set(T).intersection(Aplus))
         return node2label
     
     def roundObj(self,fracsol,basenum): #?
         """rounds objective functions(randomized way) for partition matroid base constraints
         Args:
            fracsol:   
         Returns:
            node2label:
         """   
         while True:
            node2label = {}
            for node in fracsol.keys():
                curvars = list(fracsol[node].keys())
                random.shuffle(curvars)
                p = random.random() 
                curval = 0.0
                foundvar = None
                for curvar in curvars:
                    curval += fracsol[node][curvar]
                    if curval >= p:
                       foundvar = curvar
                       break      
                if foundvar == None:
                   foundvar = curvars[-1]
                label = foundvar.replace("x{0}?".format(node),"")   
                node2label[node] = label  
            if BoundTest.checkDiffusionValidity(node2label,self.labelnodeG,"label"):
               break 
         return node2label

     def findInitialDiscreteSol(self):
         """finds initial discrete solution
         Args:
         Returns:
            initsol:
         """
         while True:
            initsol = {node: "x{0}?{1}".format(node,self.presolution[node]) for node in self.presolution.keys()} 
            for node in self.node2vars.keys():
                assert not initsol.has_key(node)
                initsol[node] = random.choice(list(self.node2vars[node]))
            if BoundTest.checkDiffusionValidity(initsol,self.labelnodeG,"var"):
               break 
         return initsol

     def delItem(self,curnode2vars,curobjval,cursol):
         """removes item from current independent set
         Args:
            curnode2vars:
            curobjval:
            cursol:
         Returns:
            bool: true/false
            cursol:
            objval:
         """
         assert len(curnode2vars.keys()) + len(self.presolution.keys()) == len(cursol.keys())
         THRES = self.EPSILON/float(math.pow(self.labelcount,4))
         foundflag = False
         for node in curnode2vars.keys():
             prelabel = cursol[node]
             if prelabel == None:
                continue 
             cursol[node] = None
             nextobjval = self.getBinaryObjVal(cursol)
             if (nextobjval-curobjval)/(-1.0*curobjval) >= THRES:
                foundflag = True
                break
             cursol[node] = prelabel
         if foundflag:    
            return True,cursol,nextobjval
         else:
            return False,cursol,curobjval
         
     def exchangeItem(self,curnode2vars,curobjval,cursol):
         """exchanges matroid items either adds or swaps to keep independent set
         Args:
            curnode2vars:
            curobjval:
            cursol:
         Returns:
            bool: true/false
            cursol:
            curobjval:
         """
         assert len(curnode2vars.keys()) + len(self.presolution.keys()) == len(cursol.keys())
         THRES = self.EPSILON/float(math.pow(self.labelcount,4))
         foundflag = False
         for node in curnode2vars.keys():
             prelabel = cursol[node]
             for var in curnode2vars[node] - set([prelabel]):
                 newlabel = var.replace("x{0}?".format(node),"") 
                 cursol[node] = newlabel
                 nextobjval = self.getBinaryObjVal(cursol)
                 if (nextobjval-curobjval)/(-1.0*curobjval) >= THRES:
                    foundflag = True
                    break
                 cursol[node] = prelabel
             if foundflag:
                break       
         if foundflag:    
            return True,cursol,nextobjval
         else:
            return False,cursol,curobjval
         
     def searchExchangeDelete(self,swapsol):
         """searches by exchange and delete operations
         Args:
            swapsol:
         Returns:
            cursol:
            curobjval:
         """
         THRES = self.EPSILON/float(math.pow(self.labelcount,4))
         newnode2vars = deepcopy(self.node2vars)
         for node in self.node2vars.keys():
             newnode2vars[node].remove(swapsol[node])
         while True:
            cursol = {node: "x{0}?{1}".format(node,self.presolution[node]) for node in self.presolution.keys()} 
            for node in newnode2vars.keys():
                assert not cursol.has_key(node)
                curvars = list(newnode2vars[node]) + [None]
                cursol[node] = random.choice(curvars) 
                curobjval = self.getBinaryObjVal(cursol)
            if curobjval != self.LOGZERO:
               break
         assert BoundTest.checkDiffusionValidity(cursol,self.labelnodeG,"var") and len(set(swapsol.keys()).symmetric_difference(set(cursol.keys()))) == 0   
         while True:
            if random.random() < 0.5:
               flag1,cursol,curobjval = self.delItem(newnode2vars,curobjval,cursol)
               if not flag1:
                  flag2,cursol,curobjval = self.exchangeItem(newnode2vars,curobjval,cursol) 
            else:
               flag1,cursol,curobjval = self.exchangeItem(newnode2vars,curobjval,cursol)
               if not flag1:
                  flag2,cursol,curobjval = self.delItem(newnode2vars,curobjval,cursol)
            if not flag2:
               break
         return cursol,curobjval 
         
     def searchOnlySwap(self):
         """searches by only swaps
         Args:
         Returns:
            cursol:
            curobjval:
         """
         THRES = self.EPSILON/float(math.pow(self.labelcount,4))
         cursol = self.findInitialDiscreteSol()
         curobjval = self.getBinaryObjVal(cursol)
         while True:
             foundflag = False
             curnodes = self.node2vars.keys()
             random.shuffle(curnodes)
             for node in curnodes:
                 prelabel = cursol[node]
                 curlabels = list(set(self.node2vars[node]).difference(set([prelabel])))
                 random.shuffle(curlabels)
                 for label in curlabels:
                     cursol[node] = label
                     nextobjval = self.getBinaryObjVal(cursol)
                     if (nextobjval-curobjval)/(-1.0*curobjval) >= THRES:
                        foundflag = True
                        break
                     cursol[node] = prelabel 
                 if foundflag:
                    break     
             if not foundflag:
                break  
             else:
                curobjval = nextobjval 
         return cursol,curobjval

     def testSolution(self,swapsol,exsol,swapobj,exobj):
         """ tests solution
         Args:
            swapsol: swap solution
            exsol: exchange solution
            swapobj: 
            exobj:
         Returns:
            bool: true/false
         """
         assert len(self.presolution.keys()) + len(self.node2vars.keys()) == len(swapsol.keys())
         for var in swapsol.values():
             assert var != None
         return True

     def findTwoBases(self,exsol):
         """finds two bases
         Args:
            exsol:
         Returns:
            base1:
            base2:
         """
         newnode2vars = deepcopy(self.node2vars)
         for node in self.newnode2vars.keys():
             if len(self.node2vars[node]) != 2:
                newnode2vars[node].remove(exsol[node]) 
         base1,base2 = {}, {}
         for node in self.newnode2vars.keys():
             posvars = list(newnode2vars[node])
             random.shuffle(posvars)
             assert len(posvars) >= 2
             base1[node] = posvars[0]
             base2[node] = posvars[1]
         return base1,base2

     def searchOptimizationTest(self,swapsol,exsol):
         """tests search optimization
         Args:
            swapsol:
            exsol:
         Returns:
            bool: true/false
         """
         assert len(swapsol.keys()) == len(self.node2vars.keys()) + len(self.presolution.keys())
         assert len(swapsol.keys()) == len(exsol.keys()) and len(set(swapsol.keys()).symmetric_difference(set(exsol.keys()))) == 0
         newnode2vars = deepcopy(self.node2vars)
         for node in self.node2vars.keys():
             newnode2vars[node].remove(swapsol[node])
         for node in newnode2vars.keys():
             assert len(newnode2vars[node]) != 0
         assert sum([len(newnode2vars[node]) for node in newnode2vars.keys()]) + len(self.node2vars.keys()) == sum([len(self.node2vars[node]) for node in self.node2vars.keys()])    
         return True
      
     def searchOptimizeObj(self):
         """search optimization based on swaps
         Args:
         Returns:
            retsol:
         """
         swapsol,objval1 = self.searchOnlySwap()
         exsol,objval2 = self.searchExchangeDelete(swapsol)
         base1,base2 = self.findTwoBases(exsol)
         newkeys = set(base1.keys()).union(set(exsol.keys()))
         for node in newkeys:
             if not base1.has_key(node):
                base1[node] = exsol[node]
             if not base2.has_key(node):
                base2[node] = exsol[node]
         assert len(base1.keys()) == len(base2.keys()) and len(base1.keys()) == len(swapsol.keys())  
         print "solinfo"
         print sol1.keys()
         print sol2.keys()
         print base1.keys()
         print base2.keys()
         exit(1)
         assert self.checkDiffusionValidity(base1) and self.checkDiffusionValidity(base2) and self.checkDiffusionValidity(swapsol)  
         sols = [swapsol,base1,base2obj]
         objs = [self.getBinaryObjVal(tsol) for tsol in sols]
         solindex = objs.index(min(objs))
         retsol,retsolobj = deepcopy(sols[solindex]), objs[solindex]
         if globals()["isTest"]:
            self.searchOptimizationTest(swapsol,exsol) 
            BoundTest.checkSolution(self.linear,self.quad,retsol,retsolobj,self.labelnodeG,"var")
            BoundTest.checkDiffusionValidity(retsol,self.labelnodeG,"var")
            labelgraph = boundInferUtil.getLabelDependency(self.node2vars,self.smodel,"normal")[0]
            boundInferUtil.checkValidityDirectly(retsol,self.G,labelgraph,"var")
         sol = {time: {state: set() for state in [Trace.SUSCEPTIBLE,Trace.EXPOSED,Trace.INFECTED,Trace.RECOVERED]} for time in xrange(self.maxtime+1)}   
         for node in retsol.keys():
             times = HistoryUtil.var2Times(retsol[node])
             for state in times.keys():
                 sol[times[state]][state].add(node) 
         return sol   
 

isTest = False    
def runner(algo,algoinfo,infermode,curstates,G,smodel,resultfile,runfolder,obstimes,TESTMODE=False):
    """infers discrete history for bounded case
       Args:
         algo: algorithm
         algoinfo: algo parameters
         infermode: infer mode
         curstates: given states information
         G: Graph
         smodel:
         resultfile: result file
         runfolder:
         obstimes: observed times of curstates
         TESTMODE:
       Returns:
         sol: 
    """
    globals()["isTest"] = TESTMODE
    assert algo in ["MatroidSub"] and algoinfo.has_key("method")
    assert checkUtil.checkPreTrace(curstates,obstimes)
    ensemble = algoinfo["ensemble"]
    enscount = 1
    if ensemble:
       enscount = algoinfo["enscount"]
    allsol = []
    for index in xrange(enscount):
        modcurstates = {obstimes[timeindex]:deepcopy(curstates[timeindex]) for timeindex in xrange(len(obstimes))}
        objclass = BoundObjFunc(G,smodel,modcurstates,algoinfo["method"])
        history = objclass.optimizeObj()
        allsol.append(deepcopy(history)) 
    if ensemble:
       history = DisEns.getEnsembleSolution(allsol,infermode,curstates,G,smodel) #??
    if infermode == "Spreader":
       foundtime = None
       for time in sorted(allsol[0].keys()):
           if len(allsol[0][time][Trace.INFECTED]) != 0:
              foundtime = time
              break 
       history = allsol[0][foundtime], 0
       print "Inferred Spreader Solution at time {0}".format(history[1])
       print history[0]
    elif infermode == "History":     
       history = allsol[0]
       print "modcurstates:"
       for time in sorted(modcurstates.keys()):  
           print time,len(modcurstates[time][Trace.SUSCEPTIBLE]),len(modcurstates[time][Trace.EXPOSED]),len(modcurstates[time][Trace.INFECTED]),len(modcurstates[time][Trace.RECOVERED])        
       print "inferred history solution"
       for time in sorted(history.keys()):
           print time,len(history[time][Trace.SUSCEPTIBLE]),len(history[time][Trace.EXPOSED]),len(history[time][Trace.INFECTED]),len(history[time][Trace.RECOVERED])
    if globals()["isTest"]:
       checkUtil.checkSolution(infermode,"bound",history,G,smodel,"detailed")   #min may not have to be at 0!
    exit(1)   
    InputOutput.writeHistoryResult(history,infermode,resultfile)   

    
if __name__ == "__main__":
    with gzip.open(sys.argv[1],"rb") as infile:
        algo = cPickle.load(infile)
        algoinfo = cPickle.load(infile)
        infermode = cPickle.load(infile)
        curstates = cPickle.load(infile)
        G = cPickle.load(infile)
        smodel = cPickle.load(infile)
        resultfile = cPickle.load(infile)
        runfolder = cPickle.load(infile)
        obstimes = cPickle.load(infile)
    runner(algo,algoinfo,infermode,curstates,G,smodel,resultfile,obstimes)
