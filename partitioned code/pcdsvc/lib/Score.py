import networkx as nx
import numpy as np
import scipy as sp
import random
import math
import myutilities as myutil
import operator
import string
import itertools

class Score():
    """ Score class
    """
    def __init__():
        return
    
    @classmethod
    def auc(self,ylist, xlist,initial,last):
        if last[0] != xlist[-1] or last[1] != ylist[-1]:
           ylist.append(last[1])
           xlist.append(last[0])  
        if initial[0] != xlist[0] or initial[1] != ylist[0]:     
           ylist.insert(0,initial[1])
           xlist.insert(0,initial[0])
        return float(np.trapz(ylist, x=xlist))   
        area = 0.0     
        sortedxlist = sorted(xlist)
        for index in range(1,len(sortedxlist)):
            prex,curx = xlist[index-1:index+1]
            prexindex,curindex = xlist.index(prex), xlist.index(curx)
            prey,cury = ylist[prexindex], ylist[curxindex]
            area += (float(curx-prex)*(cury+prey)/2.0)
        return area

    @classmethod
    def normdist(self,dist):
        mysum = sum(dist)
        return [float(elem)/mysum for elem in dist]

    @staticmethod
    def roundScore2Set(retvalues,method,methodparams,evol):
         assert method in ["epsilon"] 
         rangevalid = True
         for val in retvalues.values():
             if val > 1:
                rangevalid = False
         if not rangevalid:
            maxvalue = max(retvalues.values())
            retvalues = {key: max(retvalues[key]/float(maxvalue), 0) for key in retvalues.keys()}

         if method == "epsilon":
            epsilon = methodparams
            if evol == "static": 
               retG ={0: set([(node1,node2) for node1,node2 in retvalues.keys() if retvalues[(node1,node2)] >= 1.0 - epsilon])}
            elif evol == "dynamic":
               retG = {} 
               for node1,node2,time in retvalues.keys():
                   if retvalues[(node1,node2,time)] >= 1.0 - epsilon:
                      retG.setdefault(time,set()).add_edge(node1,node2)
            return retG
   
    @classmethod
    def kl_divergence(self,dist1,dist2):
        return sum([dist1[index]*math.log(float(dist1[index])/dist2[index],2) for index in xrange(len(dist1))])

    @classmethod
    def js_divergence(self,dist1,dist2):
        maxlen = max(len(dist1),len(dist2))
        dist1.extend([0.0] * (maxlen - len(dist1)))
        dist2.extend([0.0] * (maxlen - len(dist2)))
        half = [(dist1[index]+dist2[index])/2.0 for index in xrange(maxlen)]
        return (0.5 * kl_divergence(dist1,half)) + (0.5 * kl_divergence(dist2,half))

    @staticmethod
    def getErrorScores(truthdict,inferdict,scorenames="all"):
        return
    
    @staticmethod
    def getVectorScores(truthdict,inferdict,scorenames="all"):
        return
    
    @staticmethod
    def getBinaryScores(truthdict,inferdict,scorenames="all"):
        """
        """
        scores = {}
        tp,fp,tn,fn = 0,0,0,0
        for key in truthdict.keys():
            curtp = len(set(truthdict[key]).intersection(set(inferdict[key])))
            curfn = len(set(truthdict[key]).difference(set(inferdict[key])))
            curfp = len(set(inferdict[key]).difference(set(truthdict[key])))
            allnodes = set([node for item in truthdict[key] for node in item])
            curtn = len(allnodes) * (len(allnodes) - 1) - curfp - curtp - curfn
            tp,fp,tn,fn = tuple(map(operator.add, (tp,fp,tn,fn), (curtp,curfp,curtn,curfn)))
        scores["sen"] = float(tp)/(tp+fn)
        scores["fpr"] = float(fp)/(fp+tn)
        scores["recall"] = scores["sen"]
        scores["precision"] =0.0
        if (tp+fp) != 0:
           scores["precision"] = float(tp)/(tp+fp)
        scores["acc"] = float(tp+tn)/(tp+tn+fn+fp)
        scores["spec"] = 1.0 - scores["fpr"]
        scores["f1"], scores["f2"], scores["f01"], scores["f005"], scores["f001"]  = 0.0, 0.0, 0.0, 0.0, 0.0
        if scores["precision"] + scores["recall"] != 0.0:
           scores["f1"] = float(2.0*scores["precision"]*scores["recall"])/(scores["precision"]+scores["recall"])
           scores["f2"] = float(5.0*scores["precision"]*scores["recall"])/((4.0*scores["precision"])+scores["recall"])
           scores["f01"] =float(1.01*scores["precision"]*scores["recall"])/((0.01*scores["precision"])+scores["recall"])
           scores["f005"] =float(1.0025*scores["precision"]*scores["recall"])/((0.0025*scores["precision"])+scores["recall"])
           scores["f001"] =float(1.0001*scores["precision"]*scores["recall"])/((0.0001*scores["precision"])+scores["recall"])
        if scorenames != "all":
           for score in scores.keys(): 
               if score not in scorenames:
                  del scores[score]
        return scores
    
    @classmethod
    def getRocScore(self,xlist,ylist):
        initial,last = [0.0,0.0],[1.0,1.0]
        return auc(ylist,xlist,initial,last)

    @classmethod
    def getPrScore(self,xlist,ylist):
        initial,last = [1.0,1.0],[1.0,0.0]
        return auc(ylist,xlist,initial,last)

    @classmethod
    def getAreaScore(self,xlist,ylist,score):
        assert score in ["roc","pr","powerroc","bedroc","logroc"]
        if score == "roc":
           return getRocScore(xlist,ylist)
        elif score == "pr":
           return getPrScore(xlist,ylist)
        elif score == "powerroc":
           return getBedrocScore(xlist,ylist,"power")
        elif score == "bedroc":
           return getBedrocScore(xlist,ylist,"exponential")
        elif score == "logroc":
           return getBedrocScore(xlist,ylist,"logarithm")

    @classmethod
    def getBedrocScore(self,xlist,ylist,score):
        """bedroc is exponential transformation!
        """
        crocdir = "croc"
        codepath = "{0}/my{1}.py".format(crocdir,score)   
        tempfile = ''.join(random.choice(string.ascii_uppercase) for x in range(30))
        with open(tempfile,"w") as file:
           file.write("\n".join(["{0}\t{1}".format(xlist[index],ylist[index]) for index in xrange(len(xlist))]) + "\n")
        outfile = "out_{0}".format(tempfile)
        os.system("python {0} < {1} > {2}".format(codepath,tempfile,outfile))
        with open(outfilename,"r") as file:
           for line in file:
               score = float(line.rstrip())
        for filename in [tempfilename,outfilename]:
           os.system("rm -rf {0}".format(filename))
        return score     
