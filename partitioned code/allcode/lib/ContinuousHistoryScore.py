#Not completely done yet
#
import networkx as nx
import numpy as np
import scipy as sp
import random
import math
import sys
import os
import myutilities as myutil
import scipy.stats

class ContinuousHistoryScore():
    """ class for history infer score estimation for cont. distribution
    """
    def __init__(self):
        return

    @classmethod 
    def getPartGraphScore(self,ltruth,lsol,G,iternum,scoretype):
        """history prediction score from graph structure
        Args:
          ltruth: truth list
          lsol: our solution list
          G: Graph
          iternum: iteration count
          scoretype: scoretype
        Returns:
          score: avg graphwise closeness score 
        """
        methodname = "getPart{0}Score".format(scoretype)
        method = getattr(self,methodname)
        score = method(lsol, ltruth)
        startnodes = ltruth[0]
        for index in xrange(iternum-1):
            lnew = tracegen(startnodes,G) #??
            method = getattr(self,methodname)
            score += method(lsol, lnew)
        score /= float(len(iternum))
        return score

    @classmethod
    def getSpreaderScore(self,ltruth, lsol):
        """initial spreaders prediction score
        Args:
          ltruth: truth list
          lsol: our solution list
        Returns:
          scores: prediction scores 
        """
        scores = {}
        fn = set(ltruth[0]).difference(lsol[0])
        tp = set(ltruth[0]).intersection(lsol[0])
        fp = set(lsol[0]).difference(ltruth[0])
        tn = len(ltruth[0] + lsol[0]) - tp - fp - fn
    
        scores["sen"] = float(tp) / (tp+fn)
        scores["fpr"]= float(fp) / (fp+tn)
        scores["recall"] = sen
        scores["precision"] = 0.0
        if (tp+fp) != 0:
           scores["precision"] = float(tp) / (tp+fp)
        scores["acc"] = float(tp+tn) / (tp+tn+fn+fp)
        scores["spec"] = 1.0 - fpr #true negative rate
        scores["f1score"], scores["f2score"], scores["f01score"] = 0.0, 0.0, 0,0
        if precision+recall != 0.0:
           scores["f1score"] = float(2.0*precision*recall) / (precision+recall)
           scores["f2score"] = float(5.0*precision*recall) / ((4.0*precision)+recall)
           scores["f01score"] = float(1.01*precision*recall) / ((0.01*precision)+recall)
        return scores

    @classmethod
    def getEnsembleScore(self,ltruth, lsollist):
        """returns Ensemble score estimated by
        Args:
          ltruth: 
          lsollist: 
        Returns:
          scores: prediction scores 
        """
        scores = {}
        return scores

    @classmethod
    def getPartFootruleScore(self,l1, l2):
        """returns partial Spearman's footrule score
        Args:
           l1,l2: two lists
        Returns:
          score: partial footrule score
        """
        allnodes = reduce(lambda x,y: x.extend(y),l1)
        l1dict = {item: index for item in l1[index] for index in range(0,len(l1))}
        l2dict = {item: index for item in l2[index] for index in range(0,len(l2))}
        score = sum([abs(l1dict[node] - l2dict[node]) for node in allnodes])
        return 1.0 - score/(0.5 * len(allnodes)**2)

    @classmethod
    def getPartKenTauScore(self,l1, l2):
        """returns partial Kendall Tau score
        Args:
          l1,l2: two lists
        Returns:
          score: partial kendall tau score
        """
        l1dict = {item: index for item in l1[index] for index in range(0,len(l1))}
        l2dict = {item: index for item in l2[index] for index in range(0,len(l2))}
        allnodes = reduce(lambda x,y: x.extend(y),l1)
        concord = 0.0
        for index1 in xrange(0,len(allnodes)):
            node1 = allnodes[index1]
            for index2 in xrange(index1+1,len(allnodes)):
                node2 = allnodes[index2]
                if ((l1dict[node1] - l1dict[node2]) * (l2dict[node1] - l2dict[node2])) > 0:
                   concord += 1
                elif l1dict[node1] == l1dict[node2] and l2dict[node1] == l2dict[node2] :
                   concord += 1 
        score /= float(0.5* len(allnodes) * (len(allnodes)-1))
        return score

