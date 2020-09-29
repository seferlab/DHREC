#Plot Generator for Discrete Diffusion History Reconstruction for Unknown Case
#
import networkx as nx
import os
import sys
import random
import numpy as np
import fnmatch
import itertools
import os
sys.path.append("..")
sys.path.append("../lib")
import myutilities as myutil
import cPickle
from copy import deepcopy
from Trace import Trace
from InputOutput import InputOutput
from DiscreteTrace import DiscreteTrace as DisTrace
from ContinuousTrace import ContinuousTrace as ContTrace
from DiscreteHistoryScore import DiscreteHistoryScore as DisScore
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def getDistInfo(prob,smodel,realdata,infermode):
    """returns dists
    Args:
       prob:
       smodel:
       realdata:
       infermode:
    Returns:
       graphs:
    """                   
    prefixes = []
    if smodel == "si":
       if realdata == "real":
          prefix = "dis-si-s2i_expo_0.5_0.5_sprob_1.0_1.0"
       else:
          #prefix = "dis-si-s2i_expo_0.2_1.0_sprob_1.0_1.0" 
          #prefix = "dis-si-s2i_expo_0.2_0.2_sprob_1.0_1.0"
          prefix = "dis-si-s2i_expo_0.1_0.5_sprob_1.0_1.0"
       prefixes.append(prefix)
    elif smodel == "sir":
       if realdata == "real":
          prefix = "dis-sir-i2r_expo_0.8_1.0_s2i_weibull_9.5_2.3_sprob_0.05"
       else:
          prefix = "dis-sir-i2r_expo_0.8_1.0_s2i_expo_0.2_1.0_sprob_1.0_1.0"
       prefixes.append(prefix) 
    elif smodel == "seir":
       if realdata == "real":
          prefix = "dis-seir-i2r_expo_0.05_e2i_weibull_1.1_2.1_s2e_weibull_9.5_2.3_sprob_0.08"
       else:
          prefix = "dis-seir-i2r_expo_0.5_2.0_e2i_expo_0.5_2.0_s2e_expo_0.5_2.0_sprob_1.0_1.0"
       prefixes.append(prefix) 
    return prefixes


def retAllFiles(folder,extend):
    """returns all files under a directory
    Args:
      folder:
      extend:
    Returns:
      matches:
    """
    matches = []
    for root, dirnames, filenames in os.walk(folder):
        for filename in fnmatch.filter(filenames,extend):
            matches.append(os.path.join(root, filename))
    return matches


def selectUnderFiles(mainfolder,extend,selectmap,status):
    """ gets result paths for algo obeying selectmap
    Args:
       resfolder:
       extend:
       selectmap:
       status:
    Returns:
       retpaths:
    """
    retfiles = searchFiles(mainfolder,extend,selectmap,status)
    return retfiles


def searchFiles(mainfolder,extend,selectmap,status):
    """gets result paths for algo obeying selectmap
    Args:
       resfolder:
       extend:
       selectmap:
       status:
    Returns:
       retpaths:
    """
    assert status in ["exact","find"]
    retpaths = []
    if len(selectmap.keys()) != 0:
       for folder in myutil.listdirectories(mainfolder):
           folder2 = "{0}/{1}".format(mainfolder,folder)
           checkFolder = (status == "exact" and folder == selectmap[0]) or (status == "find" and folder.find(selectmap[0]) != -1)
           if not selectmap.has_key(0) or checkFolder:
              newmap = {index-1: selectmap[index] for index in selectmap.keys() if index != 0}
              retpaths.extend(searchFiles(folder2,extend,newmap,status))
       for filename in myutil.listfiles(mainfolder):
           checkFile = (status == "exact" and filename == selectmap[0]) or (status == "find" and filename.find(selectmap[0]) != -1)
           if selectmap.has_key(0) and len(selectmap.keys()) == 1 and checkFile and filename.endswith(extend):
              filepath = "{0}/{1}".format(mainfolder,filename)
              retpaths.append(filepath)       
    else:
       for folder in myutil.listdirectories(mainfolder):
           folder2 = "{0}/{1}".format(mainfolder,folder)
           retpaths.extend(searchFiles(folder2,extend,{},status))
       for filename in myutil.listfiles(mainfolder):
           if filename.endswith(extend):
              filepath = "{0}/{1}".format(mainfolder,filename)
              retpaths.append(filepath)
    return retpaths


def getGraphs(prob,smodel,diststr,graphfolder):
    """gets graphs
    Args:
       prob:
       smodel:
       dist:
       graphfolder:
    Returns:
       graphs:
    """
    gpaths = []
    for gfile in myutil.listfiles(graphfolder):
        gpath = "{0}/{1}".format(graphfolder,gfile)
        splitted = gfile.split("-") 
        if len(splitted) == 3:
           tprob,tsmodel,filename = splitted
           if tprob == prob and tsmodel == smodel:
              gpaths.append(gpath)
        elif len(splitted) == 4:    
           tprob,tsmodel,tdiststr,filename = splitted
           tdiststr2 = "{0}-{1}-{2}".format(tprob,tsmodel,tdiststr)
           if tprob == prob and tsmodel == smodel and tdiststr2 == diststr:
              gpaths.append(gpath)   
    return gpaths
      
def getFracs():
    """return fracs
    """ 
    return [0.15,0.26,0.34,0.4,0.51,0.6,0.76,0.9]

def getNoises():
    """return fracs
    """
    return [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

def getDistParams():
    """return dist params
    """
    return [0.1,0.2,0.25,0.5]

def getSprobParams():
    """return spread probability params
    """
    return [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

def getSpreaderCounts():
    """return initial spreader counts
    """
    return [1,2,3,4,5,6,7,8,9,10]

def getSnapshots():
    """return sample count(number of samples)
    """
    return [1,2,3,4,5,6,7,8]

def getPlotnameLabels(infermode): 
    """return plotname labels
    Args:
       infermode:
    """   
    labelmap = {"sprobparam":"Spreading Probability","distparam":"Distribution Parameters","snapshots":"# Snapshots ("+r'$|T_{D}|$'+")","spreadercount":"# Initial Spreaders", "noise": "Noise Rate"}
    if infermode == "Spreader":
       labelmap["fracs"] = "Min Snapshot Ratio ("+r'$f_{\min}$'+")"
    elif infermode == "History":
       labelmap["fracs"] = "Max Snapshot Ratio ("+r'$f_{\max}$'+")"
    return labelmap
        
def getScoremapLabels():
    """ return scoremap labels
    """    
    return {"F1":"F1", "F01":"F0.1", "acc":"Accuracy", "KenTau": r'1 - $\tau$', "KenTauCorrelation": r'$\tau_{B}$',"Hausdorff":"1 - H", "Footrule":"1 - F", ("Graph","jaccardsim"): r'$J_{G}$', ("Graph","match"): r'$\overline{M}_{G}$', ("Graph","rand"): r'$\mathcal{R}_{G}$', "jaccardsim": 'J', "jaccarddis":r'$J_{\delta}$', "rand":r'$\mathcal{R}$'}    

def getScoreLimits():
    """
    """   
    return {"F1":(0,1.0),"F01":(0,1.0),"acc":(0,1.0),"KenTau":(0,1.0),"Footrule":(0,1.0),"KenTauCorrelation":(-1.0,1.0),"Hausdorff":(0,1.0),"jaccardsim":(0,1.0),"jaccarddis":(0,1.0),"rand":(0,1.0),"match":(0,1.0)}

def getScore(truth,sol,graphparams,scoreparams,obstimes):    
    """returns score
    Args:
      truth:
      sol:
      graphparams:
      scoreparams:
      obstimes:
    Returns:
      score: 
    """
    if min(truth.keys()) < 0:
       assert max(obstimes) == -1 * min(truth.keys())
    G,evol,smodel,infermode,inter = graphparams
    curscore,scoreparam,usestates = scoreparams
    if type(curscore) == tuple and curscore[0] == "Graph":
       if len(sol.values()[0][Trace.INFECTED]) == 0:
          score = 0.0
       elif curscore[1] == "match": 
          score = DisScore.getMatchingScore(truth,sol,G,usestates)
          assert score >= 0.0 and score <= 1.0
       else:    
          score = DisScore.getPartGraphScore(truth,sol,G,10,curscore[1],smodel,evol,obstimes,inter)
    if curscore == "KenTau":
       scoremethod = "get{0}Score".format(curscore) 
       score =  getattr(DisScore,scoremethod)(truth,sol,0.5,usestates)
    elif curscore in ["Footrule","Hausdorff","KenTauCorrelation"]: 
       scoremethod = "get{0}Score".format(curscore)
       score = getattr(DisScore,scoremethod)(truth,sol,usestates)  
    elif curscore in ["sen","fpr","recall","precision","acc","spec","F1","F01","F01","jaccardsim","jaccarddis"]:
       scoremethod = "getPredictionScore"    
       score = getattr(DisScore,scoremethod)(truth,sol)[curscore]
    elif curscore == "roc":
       scoremethod = "getPredictionScore"    
       rocscores = getattr(DisScore,scoremethod)(truth,sol)
       return (rocscores["fpr"],rocscores["sen"])
    return score

def getTraceMap(plotname,fixedname,noisetype,tracefile,srate,x,traceside): 
    """
    Args:

    Returns:
       tracemap:
    """
    tracemap = {} 
    if plotname in ["fracs", "snapshots"]:
       noise,startcount = traceside
       tracemap = {0:fixedname, 1:"{0}-{1}-{2}".format(noisetype,noise,srate), 2:str(startcount)}
    elif plotname == "spreadercount":
       noise = traceside[0]
       tracemap = {0:fixedname, 1:"{0}-{1}-{2}".format(noisetype,noise,srate), 2:str(x)}
    elif plotname == "noise":
       startcount = traceside[0] 
       tracemap = {0:fixedname, 1:"{0}-{1}-{2}".format(noisetype,x,srate), 2:str(startcount)}
    elif plotname == "distparam":
       pass
    elif plotname == "sprobparam":
       pass
    if tracefile != None:
       tracemap[3] = tracefile
    return tracemap   

def getResultMap(plotname,fixedname,noisetype,curtfile,infermode,srate,x,resultside):
    """
    Args:

    Returns:
    """
    if plotname == "fracs":
       noise,startcount,timecount = resultside
       fracpoints = [(x*index)/timecount for index in xrange(1,timecount+1)]
       resmap = {0:fixedname, 1:"{0}-{1}-{2}-{3}-{4}".format(noisetype,noise,srate,startcount,curtfile), 2:"frac{0}-count{1}".format(x,timecount), 3:infermode}
    elif plotname == "snapshots":
       noise,startcount,timefrac = resultside 
       fracpoints = [(timefrac*index)/x for index in xrange(1,x+1)]
       resmap = {0:fixedname, 1:"{0}-{1}-{2}-{3}-{4}".format(noisetype,noise,srate,startcount,curtfile), 2:"frac{0}-count{1}".format(timefrac,x), 3:infermode}   
    elif plotname == "spreadercount":
       noise,timecount,timefrac = resultside 
       fracpoints = [(timefrac*index)/timecount for index in xrange(1,timecount+1)]
       resmap = {0:fixedname, 1:"{0}-{1}-{2}-{3}-{4}".format(noisetype,noise,srate,x,curtfile), 2:"frac{0}-count{1}".format(timefrac,timecount), 3:infermode}
    elif plotname == "noise":
       startcount,timecount,timefrac = resultside 
       fracpoints = [(timefrac*index)/timecount for index in xrange(1,timecount+1)]
       resmap = {0:fixedname, 1:"{0}-{1}-{2}-{3}-{4}".format(noisetype,x,srate,startcount,curtfile), 2:"frac{0}-count{1}".format(timefrac,timecount), 3:infermode}
    elif plotname == "distparam":
       pass
    elif plotname == "sprobparam":
      pass
    return resmap,fracpoints
      


def prepareData(algoparams,traceparams,plotparams,sideparams):
    """prepares data for plotting
    Args:
      algoparams:
      traceparams:
      plotparams:
      sideparams:
    Returns:
      valmap:
    """
    [algos,G,fixedname,infermode,inter,algomap] = algoparams
    [curscore,plotname,usestates,scoreparam,plotpath] = plotparams
    [tracefolder,resultfolder,srate,noisetype,evol,smodel,prob,tracefile] = traceparams
    assert plotname in ["sprobparam", "distparam", "fracs", "snapshots", "spreadercount", "noise"]
    if plotname == "fracs":
       [xs,noise,timecount,startcount] = sideparams
       traceside = [noise,startcount]
       resultside = [noise,startcount,timecount]
    elif plotname == "snapshots":
       [xs,timefrac,noise,startcount] = sideparams
       traceside = [noise,startcount]
       resultside = [noise,startcount,timefrac]
    elif plotname == "spreadercount":
       [xs,noise,timefrac,timecount] = sideparams
       traceside = [noise]
       resultside = [noise,timecount,timefrac]
    elif plotname == "noise":
       [xs,timefrac,timecount,startcount] = sideparams
       traceside = [startcount]
       resultside = [startcount,timecount,timefrac]
    elif plotname == "sprobparam":
       pass
    elif plotname == "distparam":
       pass
    
    realalgos = set()
    for realalgo in set(algomap.values()):
        flag = False
        for algo in algos:
            if algomap[algo] == realalgo:
               flag = True
               break 
        if flag:
           realalgos.add(realalgo) 
    valmap = {algo: {"x":[], "y":[] } for algo in realalgos}
    opmap = {"KenTau":max, "KenTauCorrelation":max, "Footrule":max, "Hausdorff":max, "F1":max, "F01":max, "acc":max, "match":max, "jaccardsim":max, "jaccarddis":min, "roc":max, "rand":max}
    if curscore != "KenTauCorrelation": 
       op2add = {min:1.0,max:0.0}
    else:
       op2add = {min:1.0,max:-1.0}    
    algo2y = {algo: {x:[] for x in xs} for algo in realalgos}
    if curscore == "roc":
       algo2x = {algo: {x:[] for x in xs} for algo in realalgos}
    if type(curscore) == tuple:
       op = opmap[curscore[1]]
    else:
       op = opmap[curscore]
    for x in xs:
        extensionstr = fixedname
        tracemap = getTraceMap(plotname,fixedname,noisetype,tracefile,srate,x,traceside)
        tracepaths = selectUnderFiles(tracefolder,'.plain',tracemap,"exact")
        print tracefolder
        print tracemap
        print x,len(tracepaths)
        for tracepath in tracepaths:
            curtfile = tracepath.split("/")[-1]
            trace = InputOutput.readPlainTrace(tracepath,prob)
            resmap,fracpoints = getResultMap(plotname,fixedname,noisetype,curtfile,infermode,srate,x,resultside) 
            giventimes = list(set([int(round(frac * Trace.getMaxTraceTime(trace))) for frac in fracpoints]))
            #print "giventimes",giventimes
            truth = DisScore.trace2History(trace,evol,smodel,giventimes)
            assert -1 * min(truth.keys()) == max(giventimes) and max(truth.keys()) == 0
            temptruth = deepcopy(truth)
            zeroref = -1*min(truth.keys())
            truth = {prekey+zeroref: temptruth[prekey] for prekey in temptruth.keys()}
            for algo in algos:
                resmap[4] = algo
                #print "info"
                #print algo
                #print resultfolder
                #print resmap
                resultfiles = selectUnderFiles(resultfolder,'.hist',resmap,"exact")
                #print "info"
                #print resultfolder
                #print len(resultfiles)
                #print resmap
                #exit(1)
                for curresultpath in resultfiles:
                    print curresultpath
                    if infermode == "Spreader":
                       mintime,sol = InputOutput.readHistoryResult(curresultpath,infermode)
                       sol = {mintime: sol}
                    elif infermode == "History":
                       sol = InputOutput.readHistoryResult(curresultpath,infermode)
                    if inter == "bound":
                       if infermode == "Spreader":
                          assert min(sol.keys()) == min(truth.keys()) 
                       elif infermode == "History":  
                          assert len(set(sol.keys()).symmetric_difference(truth.keys())) == 0
                    for scorestate in usestates:
                        for time in truth.keys():
                            if not truth[time].has_key(scorestate):
                               truth[time][scorestate] = set() 
                        for time in sol.keys():
                            if not sol[time].has_key(scorestate):
                               sol[time][scorestate] = set()  
                    graphparams = [G,evol,smodel,infermode,inter]
                    scoreparams = [curscore,scoreparam,usestates]
                    if curscore != "roc":
                       tscore = getScore(truth,sol,graphparams,scoreparams,giventimes) 
                       if curscore in ["Hausdorff","KenTau","Footrule"]: 
                          algo2y[algomap[algo]][x].append(1.0-tscore)
                       else:
                          algo2y[algomap[algo]][x].append(tscore)    
                    else:
                       curfpr,cursen = getScore(truth,sol,graphparams,scoreparams)
                       algo2x[algomap[algo]][x].append(curfpr) 
                       algo2y[algomap[algo]][x].append(cursen)
    if curscore != "roc":
       for realalgo in algo2y.keys():
           for x in sorted(algo2y[realalgo].keys()):
               valmap[realalgo]["x"].append(x)
               valmap[realalgo]["y"].append(op(algo2y[realalgo][x] + [op2add[op]]))
    else:
       pass 
    return valmap

def assignMarker(valmap):
    """assigns markers
    Args:
      valmap:
    Returns:
      symmap: 
    """    
    #colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
    #symbols
    #'^' : triangle up
    #'>' : triangle right
    #'v' : triangle down
    #'<' : triangle left
    #'d' : diamond
    #'p' : pentagram
    #'h' : hexagon
    #'8' : octagon
    globmarkers = {"NetSleuth":"d","DHR-sub":"s","DHR-sub-ens":"p","DHR-pcvc":"*","DHR-pcvc-ens":"x","DHR-pcdsvc":"o","DHR-pcdsvc-ens":"h","Rumor":"+", "KEffectors":".","GreedyForward":"^"}
    globcolors = {"NetSleuth":"b","DHR-sub":"k","DHR-sub-ens":"m","DHR-pcvc":"c","DHR-pcvc-ens":"c","DHR-pcdsvc":"y","DHR-pcdsvc-ens":"y","Rumor":"g", "KEffectors":"r","GreedyForward":"m"}
    globlabels = {"NetSleuth":"NetSleuth","DHR-sub":"DHR-sub","DHR-sub-ens":"DHR-sub-ens","DHR-pcvc":"DHR-pcvc","DHR-pcvc-ens":"DHR-pcvc-ens","DHR-pcdsvc":"DHR-pcdsvc","DHR-pcdsvc-ens":"DHR-pcdsvc-ens","Rumor":"Rumor","KEffectors":"KEffectors","GreedyForward":"GreedyForward"}
    symmap = {"marker":{}, "colors":{}, "labels":{}}
    inmap,index = {0:"marker", 1:"colors", 2:"labels"}, 0   
    for cursym in [globmarkers,globcolors,globlabels]:
        for algo in valmap.keys():
            symmap[inmap[index]][algo] = cursym[algo]
        index += 1
    return symmap

def plotHeatMap(data,indexmap,curscore,plotpath,axlabels,infermode):
    """plots heatmap and saves it
    Args: 
       data: data in array format
       indexmap: script mapper
       curscore: score
       plotpath: plotpath
       axlabels: axis labels
       infermode:
    Returns:
    """
    if infermode == "History":
       labelmap = {"frac":"Max Snapshot Ratio ("+r'$f_{\max}$'+")","timecount":"# Snapshots ("+r'$|T_{D}|$'+")","startcount":"# Initial Spreaders"}
    elif infermode == "Spreader":
       labelmap = {"frac":"Min Snapshot Ratio ("+r'$f_{\min}$'+")","timecount":"# Snapshots ("+r'$|T_{D}|$'+")","startcount":"# Initial Spreaders"} 
    plt.clf()
    plt.rc('font', family='serif', size=25)
    FSIZE = 24
    DPI = 1000
    fig, ax = plt.subplots()
    fig.set_size_inches(10,8)
    mycmap = plt.cm.Reds #plt.cm.Blues
    heatmsap = ax.pcolor(data,cmap=mycmap,alpha=0.8,edgecolors='k',linewidths=1)
    ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)   
    ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)
    xlabels = [indexmap["row"][index] for index in sorted(indexmap["row"].keys())]
    ylabels = [indexmap["col"][index] for index in sorted(indexmap["col"].keys())]
    #print xlabels
    #print ylabels
    #print np.shape(data)
    ax.set_xticklabels(ylabels, minor=False, fontsize=FSIZE)
    ax.set_yticklabels(xlabels, minor=False, fontsize=FSIZE)
    ax.grid(False)
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks(): 
       t.tick1On = False 
       t.tick2On = False 
    for t in ax.yaxis.get_major_ticks(): 
       t.tick1On = False 
       t.tick2On = False
    #minscore,maxscore = getScoreLimits()[curscore] 
    #barticks = [minscore + (index * (maxscore-minscore)/2) for index in xrange(3)]
    #barstrs = [str(tick) for tick in barticks]
    #cbar = fig.colorbar(heatmsap, ticks=barticks)
    #cbar.ax.set_yticklabels(barstrs)
    #fig.set_size_inches(9,9)
    cb = fig.colorbar(heatmsap)
    plt.text(0.5,1.05,getScoremapLabels()[curscore],horizontalalignment='center',fontsize=FSIZE*2,transform = ax.transAxes)
    #ax.set_title(getScoremapLabels()[curscore],fontsize=FSIZE)
    ax.set_xlabel(labelmap[axlabels[0]],fontsize=FSIZE)
    ax.set_ylabel(labelmap[axlabels[1]],fontsize=FSIZE)
    plt.savefig(plotpath, dpi=DPI)
    
    
def genHeatMap(algoparams,traceparams,plotparams,sideparams):
    """ generates heatmap
    Args:
      algoparams:
      traceparams:
      plotparams:
      sideparams:
    Returns:
    """
    [algos,G,fixedname,infermode,inter,algomap] = algoparams
    [curscore,plotname,usestates,scoreparam,plotpath] = plotparams
    [tracefolder,resultfolder,srate,noisetype,evol,smodel,prob,tracefile] = traceparams
    assert plotname in ["timecount-frac","frac-timecount","frac-startcount","timecount-startcount"]
    assert len(set(algomap[algo] for algo in algos)) == 1
    realalgo = list(set(algomap[algo] for algo in algos))[0]
    sentalgomap = {}
    for algo in algomap.keys():
        if algomap[algo] == realalgo:
           sentalgomap[algo] = algomap[algo]           
    myalgoparams = [algos,G,fixedname,infermode,inter,sentalgomap]
    if plotname == "timecount-frac":
       timecounts,fracs,noise,startcount = sideparams
       valmap = {}
       for frac in fracs:
           mysideparams = [timecounts,frac,noise,startcount]
           myplotparams = [curscore,"snapshots",usestates,scoreparam,plotpath]
           valmap[frac] = prepareData(myalgoparams,traceparams,myplotparams,mysideparams)
       collist, rowlist = sideparams[0:2]  
    elif plotname == "frac-timecount":
       fracs,timecounts,noise,startcount = sideparams
       valmap = {}
       for timecount in timecounts:
           mysideparams = [fracs,noise,timecount,startcount]
           myplotparams = [curscore,"fracs",usestates,scoreparam,plotpath]
           valmap[timecount] = prepareData(myalgoparams,traceparams,myplotparams,mysideparams)
       collist, rowlist = sideparams[0:2]
    elif plotname == "frac-startcount":
       fracs,startcounts,noise,timecount = sideparams
       valmap = {}
       for startcount in startcounts:
           mysideparams = [fracs,noise,timecount,startcount]
           myplotparams = [curscore,"fracs",usestates,scoreparam,plotpath]
           valmap[startcount] = prepareData(myalgoparams,traceparams,myplotparams,mysideparams)
       collist, rowlist = sideparams[0:2]
    elif plotname == "timecount-startcount":
       timecounts,startcounts,noise,timefrac = sideparams
       valmap = {}
       for startcount in startcounts:
           mysideparams = [timecounts,timefrac,noise,startcount]
           myplotparams = [curscore,"snapshots",usestates,scoreparam,plotpath]
           valmap[startcount] = prepareData(myalgoparams,traceparams,myplotparams,mysideparams)
       collist, rowlist = sideparams[0:2]          
    assert sorted(collist) == collist and sorted(rowlist) == rowlist   
    print "info ",plotpath,plotname   
    print collist
    print rowlist
    singlealgo = list(set(algo for timecount in valmap.keys() for algo in valmap[timecount].keys()))[0]
    for timecount in sorted(valmap.keys()):
        for index in xrange(len(collist)):
            col = collist[index]
            print timecount,col,valmap[timecount][singlealgo]["y"][index]
    indexmap = {"col": {index: collist[index]  for index in xrange(len(collist))}, "row": {index: rowlist[index] for index in xrange(len(rowlist))}}   
    data = np.zeros((len(rowlist),len(collist)),dtype=np.float)
    for row in xrange(len(rowlist)):
        rowval = rowlist[row]
        for col in xrange(len(collist)):
            data[row,col] = valmap[rowval][realalgo]["y"][col]  
    axlabels = plotname.split("-")  
    plotHeatMap(data,indexmap,curscore,plotpath,axlabels,infermode)

    
def genPlot(algoparams,traceparams,plotparams,sideparams):
    """generates classical plots
    Args:
      algoparams:
      traceparams:
      plotparams:
      sideparams:
    Returns:
    """
    [algos,G,fixedname,infermode,inter,algomap] = algoparams
    [curscore,plotname,usestates,scoreparam,plotpath] = plotparams
    [tracefolder,resultfolder,srate,noisetype,evol,smodel,prob,tracefile] = traceparams
    assert plotname in ["sprobparam", "distparam", "fracs", "snapshots", "spreadercount", "noise"] 
    if plotname == "fracs":
       [xs,noise,timecount,startcount] = sideparams
    elif plotname in ["snapshots"]:
       [xs,timefrac,noise,startcount] = sideparams
    elif plotname == "spreadercount":
       [xs,noise,timefrac,timecount] = sideparams
    elif plotname == "noise":
       [xs,timefrac,timecount,startcount] = sideparams
    elif plotname == "sprobparam":
       pass
    elif plotname == "distparam":
       pass
    valmap = prepareData(algoparams,traceparams,plotparams,sideparams)
    symbolmap = assignMarker(valmap)
    
    plt.clf()
    plt.rc('font', family='serif', size=25)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13,10)
    FSIZE = 30
    YFSIZE = 60
    LEGENDSIZE = 30 #20
    MARKERSIZE = 15.0 #10.0
    DPI = 1000
    ylimit = getScoreLimits()
    if type(curscore) == tuple: 
       plt.ylim(ylimit[curscore[1]])
       ylimscore = ylimit[curscore[1]]
    else:
       plt.ylim(ylimit[curscore])
       ylimscore = ylimit[curscore]    
    xlabels = getPlotnameLabels(infermode)
    scoremap = getScoremapLabels()
    plotparts = plotpath.split("/")
    scorepath = "/".join(plotparts[0:-1]) + "/" + plotparts[-1].replace(".png",".txt")
    avgsum = sum([val for algo in valmap.keys() for val in valmap[algo]["y"]])   
    if avgsum <= sum(ylimscore) / 2.0:
       locpos = 1 
    else:    
       locpos = 4
                 
    if curscore != "roc":
       plt.xlabel(xlabels[plotname],fontsize=FSIZE)
       plt.ylabel(scoremap[curscore],fontsize=YFSIZE)
    else:
       plt.xlabel("FPR",fontsize=FSIZE)
       plt.ylabel("Sensitivity",fontsize=YFSIZE)
    with open(scorepath,"w") as outfile:
        if curscore != "roc":
           allx = set(x for algo in valmap.keys() for x in valmap[algo]["x"]) 
           xlimdict = {"fracs":0.1, "sprobparam":0.1, "distparam":1.0, "snapshots":1, "spreadercount":1, "noise":0.1}
           plt.xlim([0,max(list(allx)) + xlimdict[plotname]])
        else:
           pass
        print valmap.keys()
        for algo in valmap.keys():
            print "algo info"
            print algo
            print valmap[algo]["x"]
            print valmap[algo]["y"]
            plt.plot(valmap[algo]["x"],valmap[algo]["y"],marker=symbolmap["marker"][algo],markersize=MARKERSIZE,linestyle='None',color=symbolmap["colors"][algo],label=symbolmap["labels"][algo])
            xstr = " ".join([str(item) for item in valmap[algo]["x"]])
            ystr = " ".join([str(item) for item in valmap[algo]["y"]])
            outfile.write("{0}\t{1}\t{2}\n".format(algo,xstr,ystr))
    plt.legend(loc=locpos,prop={'size':LEGENDSIZE})
    plt.savefig(plotpath, dpi=DPI)


def getAlgos(prob,inter,infermode,smodel):
    if prob == "dis" and inter == None:
       if infermode == "Spreader":
          algoblocks = [("GreedySubCoverSingle",{"iter":None,"objmethod":"log","ensemble":False}),("GreedySubCoverSingle",{"iter":None, "objmethod":"log","ensemble":True,"enscount":2}),("GreedySubCoverSingle",{"iter":None, "objmethod":"log","ensemble":True,"enscount":3})]
          #algoblocks = [("FracCover",{"iter":None, "objmethod":"log","ensemble":False,"roundmethod":"random"})]
          #algoblocks = []
          if smodel in ["si","sir"]:
             appralgos = [("Pcdsvc",{"ensemble":False}),("Pcdsvc",{"ensemble":True,"enscount":2}),("Pcdsvc",{"ensemble":True,"enscount":5})]
             appralgos2 = [("Pcvc",{"ensemble":False}),("Pcvc",{"ensemble":True,"enscount":2}),("Pcvc",{"ensemble":True,"enscount":5})]
             algoblocks.extend(appralgos)
             algoblocks.extend(appralgos2)
          elif smodel == "seir":
             appralgos = [("MinCut",{"ensemble":False}),("MinCut",{"ensemble":True,"enscount":2}),("MinCut",{"ensemble":True,"enscount":5})]  
             algoblocks.extend(appralgos)
          #algoblocks = []
          #otheralgos = [("NetSleuth",{}),("RumorCentrality",{}),("KEffectors",{})]
          #otheralgos = [("NetSleuth",{})]
          otheralgos = [("NetSleuth",{}),("RumorCentrality",{})]
          algoblocks.extend(otheralgos)
       elif infermode == "History":
          algoblocks = [("GreedySubCoverSingle",{"iter":None, "objmethod":"log","ensemble":False})]
          #algoblocks = [("FracCover",{"iter":None, "objmethod":"log","ensemble":False,"roundmethod":"random"})]
          #algoblocks = []
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

def getAlgoMap():
    """returns algo map
    Args:
    Returns:
       algomap:   
    """  
    #algomap["FracCover-random-log-False-None"] = "DHR-con"        
    algomap = {}
    algomap["NetSleuth"] = "NetSleuth"
    algomap["RumorCentrality"] = "Rumor"
    algomap["KEffectors"] = "KEffectors"
    algomap["GreedySubCoverSingle-log-False-None"] = "DHR-sub"
    #algomap["GreedySubCoverSingle-log-True-None-5"] = "DHR-sub-ens"
    algomap["GreedySubCoverSingle-log-True-None-2"] = "DHR-sub-ens"
    algomap["GreedySubCoverSingle-log-True-None-3"] = "DHR-sub-ens"
    algomap["MinCut-False"] = "DHR-pcvc"
    algomap["MinCut-True-5"] = "DHR-pcvc-ens"
    algomap["MinCut-True-2"] = "DHR-pcvc-ens"
    algomap["MinCut-True-3"] = "DHR-pcvc-ens"
    algomap["Pcdsvc-False"] = "DHR-pcdsvc" 
    algomap["Pcdsvc-True-2"] = "DHR-pcdsvc-ens" 
    algomap["Pcdsvc-True-3"] = "DHR-pcdsvc-ens" 
    algomap["Pcdsvc-True-5"] = "DHR-pcdsvc-ens" 
    algomap["Pcvc-False"] = "DHR-pcvc" 
    algomap["Pcvc-True-2"] = "DHR-pcvc-ens"
    algomap["Pcvc-True-3"] = "DHR-pcvc-ens" 
    algomap["Pcvc-True-5"] = "DHR-pcvc-ens"
    algomap["GreedyForward-False"] = "GreedyForward" 
    #algomap["Qsap-reliability-False"] = "DHRi-min"
    #algomap["Qsap-random-ratio-False"] = "DHRi-max"
    #algomap["MatroidSub-False"] = "DHRi-sub"
    #algomap["MatroidSub-True-5"] = "DHRi-sub-ens"
    #algomap["Independent-False"] = "Independent"
    return algomap


def getAlgoStrs(plottype,prob,inter,infermode,smodel):
    """
    """
    if inter == None:
       if plottype in ["sprobparam","distparam","fracs","snapshots","spreadercount","noise"]:   
          algos = [returnAlgoStr(algo[0],algo[1]) for algo in getAlgos(prob,inter,infermode,smodel)]
       elif plottype in ["timecount-frac","frac-timecount","frac-startcount","timecount-startcount","heat-param"]:
          if infermode == "History": 
             algos = ["GreedySubCoverSingle-log-False-None"]
          elif infermode == "Spreader": 
             algos = ["GreedySubCoverSingle-log-False-None"] # 'NetSleuth', 'Rumor', 'KEffectors']     
    elif inter == "bound":
       if plottype in ["sprobparam","distparam","fracs","snapshots","spreadercount","noise"]:   
          algos = [returnAlgoStr(algo[0],algo[1]) for algo in getAlgos(prob,inter,infermode,smodel)]
       elif plottype in ["frac-timecount","frac-startcount","timecount-startcount","heat-param"]:
          algos = ["Qsap-random-ratio-False"] #"Qsap-reliability-False"
    return algos


def retAlgos(plottype,prob,inter,infermode,smodel):
    """returns algos
    Args:
       plottype:
       prob:
       inter:
       infermode:
       smodel:
    Returns:
       algos:
    """
    algos = getAlgoStrs(plottype,prob,inter,infermode,smodel)  
    return algos


def makeCheck(infermode,curscore,plottype,algos,algomap,usestates):
    """
    """    
    if infermode == "Spreader":
       assert curscore not in ["KenTau","Footrule","Hausdorff","KenTauCorrelation"]
       assert len(usestates) == 1 and usestates[0] == Trace.INFECTED 
    elif infermode == "History":
       assert curscore not in ["Graph","F1","F01","jaccardsim","jaccarddis","roc"]
    if type(curscore) == tuple:
       assert curscore[0] == "Graph" and curscore[1] in ["match","jaccardsim","jaccarddis","rand"]    
    #if curscore == "roc":
    #   assert plottype in "snapshots"
    #print algos
    #print algomap.keys()
    #print set(algos).difference(set(algomap.keys()))
    assert len(set(algos).difference(set(algomap.keys()))) == 0

    
# Plots
# 1- Spreader:
#    y: score, x: frac of max point 
#    y: score, x: number of snapshots(only when diffusion length is known)
#    y: score  x: parameter of the distribution 
#    y: score  x: initial spreader count
#    y: score  x: noise degree
# 2- History: 
#    y: score, x: frac of max point 
#    y: score, x: number of snapshots 
#    y: score  x: parameter of the distribution
#    y: score  x: noise degree
#
# roc curve (inferrability) vs number of snapshots(not exactly roc curve)
# heatmap plot dmc parameters or dist parameters(2 dimension)
# heatmap plot fracs and number of snapshots (combination of two plots above)
# runtime table: 

def main():
    """main method for plotting
    """
    plotpref = "plot"
    graphpref = "graphs"
    resultpref = "result"
    tracepref = "traces"
    realdata = "syn"
    noise = 0.0
    srate = 0
    startcount = 1
    timefrac = 0.26 #[0.15,0.26,0.34,0.51,0.6,0.76,0.9 #percentage of diffusion maxtime(peak time gibi)
    timecount = 1 #3,5,10 
    noisetype = "StateChange"
    evol = "static"
    smodel = "si"
    usestates = [Trace.INFECTED] #states to score estimated
    tracefile = "0.plain" #"0.plain" #None #"1.plain" #None
    inter = None #None
    infermode = "History" #"History" #"History"
    plottype = "fracs" #"spreadercount" #"timecount-frac" #"timecount-startcount" #"timecount-frac" #"snapshots" #"timecount-frac" #"frac-timecount" #"spreadercount" #"timecount-startcount" #"spreadercount" #"timecount-startcount" #"distparam", "fracs", "snapshots", "spreadercount", "noise", "frac-timecount", "frac-startcount", 
    resultfolderpref = ".." #"../store/storeresults"
    graphfolderpref = ".."
    tracefolderpref = ".."
    plotfolderpref = "."
    prob = "dis" #"cont"
    if infermode == "Spreader":
       #curscores = ["jaccardsim","jaccarddis","acc","F1","F01"], ("Graph", "rand"),("Graph", "jaccardsim")  "roc" curscores = [("Graph", "rand"),"acc"] 
       curscores = [("Graph", "match")] 
    elif infermode == "History":
       curscores = ["KenTauCorrelation"] #curscores = ["KenTauCorrelation","Hausdorff","KenTau","Footrule"]
    scoreparam = []
    usestatesstr = "_".join(usestates)
       
    algomap = getAlgoMap()
    algos = retAlgos(plottype,prob,inter,infermode,smodel)
    [makeCheck(infermode,curscore,plottype,algos,algomap,usestates) for curscore in curscores]
       
    graphfolder = "{0}/{1}_{2}_{3}".format(graphfolderpref,realdata,evol,graphpref)
    resultfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(resultfolderpref,resultpref,realdata,evol,smodel,prob,str(inter))
    tracefolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}".format(tracefolderpref,tracepref,realdata,evol,smodel,"edge",prob)
    plotfolder = "{0}/{1}_{2}_{3}_{4}_{5}_{6}/{7}".format(plotfolderpref,plotpref,realdata,evol,smodel,prob,str(inter),infermode)
    if not os.path.exists(plotfolder):
       os.makedirs(plotfolder)  
    filepaths = [filename for diststr in getDistInfo(prob,smodel,realdata,infermode) for filename in getGraphs(prob,smodel,diststr,graphfolder) if filename.find("deg10_")==-1]
 
    #filepaths = [filepath.replace(".edgelist",".pickle")for filepath in filepaths]
    
    #filepaths = [filepath for filepath in filepaths if filepath.find("_undirected_ff")!=-1]
    if realdata == "real":
       filepaths = [filepath for filepath in filepaths if filepath.find("sn1.pickle")!=-1]
    for filepath,curscore in list(itertools.product(filepaths,curscores)):
        #if not filepath.split("_")[-1].startswith("1."):
        #   continue 
        fixedname = filepath.split("/")[2]
        G = InputOutput.readGraphAndParams(filepath)
        curplotfolder = "{0}/{1}".format(plotfolder,fixedname) 
        if not os.path.exists(curplotfolder):
           os.makedirs(curplotfolder)
        algoparams = [algos,G,fixedname,infermode,inter,algomap]
        traceparams = [tracefolder,resultfolder,srate,noisetype,evol,smodel,prob,tracefile]
        plotparams = [curscore,plottype,usestates,scoreparam]
        
        if plottype == "fracs":
           xs = sorted(getFracs())
           if infermode == "Spreader":
              #xs = [0.15, 0.26, 0.34, 0.4, 0.51, 0.6, 0.76, 0.9]
              xs = [0.15,0.26,0.34,0.45,0.55,0.65,0.76,0.85,0.9,0.95]
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.png".format(infermode,curscore,smodel,tracefile,startcount,timecount,noise,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = [xs,noise,timecount,startcount]
        elif plottype == "noise":
           xs = sorted(getNoises()) 
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.png".format(infermode,curscore,smodel,tracefile,startcount,timecount,timefrac,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = [xs,timefrac,timecount,startcount]
        elif plottype == "spreadercount":
           xs = sorted(getSpreaderCounts())
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.png".format(infermode,curscore,smodel,tracefile,timecount,timefrac,noise,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = [xs,noise,timefrac,timecount]
        elif plottype == "snapshots":
           xs = sorted(getSnapshots())
           xs = [1,2,3,4,5,6]
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.png".format(infermode,curscore,smodel,tracefile,startcount,noise,timefrac,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = [xs,timefrac,noise,startcount]
        elif plottype == "distparam":
           xs = sorted(getDistParams())
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}.png".format(infermode,curscore,smodel,tracefile,startcount,timecount,noise,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = []
        elif plottype == "sprobparam":
           xs = sorted(getSprobParams())
           plotfile = ""
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = []  
        elif plottype == "timecount-frac": 
           xs = sorted(getSnapshots())
           ys = sorted(getFracs())
           xs = [1,2,3,5,7]
           ys = [0.26,0.34,0.51,0.76,0.9] 
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}.png".format(infermode,curscore,smodel,tracefile,startcount,noise,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = [xs,ys,noise,startcount]       
        elif plottype == "frac-timecount": 
           xs = sorted(getFracs())
           ys = sorted(getSnapshots())
           xs = [0.34,0.6,0.76,0.8,0.9]
           ys = [1,2,3,5,7]
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}.png".format(infermode,curscore,smodel,tracefile,startcount,noise,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = [xs,ys,noise,startcount]
        elif plottype == "frac-startcount":
           xs = sorted(getFracs())
           ys = sorted(getSpreaderCounts())
           xs = [0.15,0.26,0.34,0.51,0.6]
           ys = [1,2,3,5,7]
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}.png".format(infermode,curscore,smodel,tracefile,timecount,noise,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = [xs,ys,noise,timecount]
        elif plottype == "timecount-startcount":
           xs = sorted(getSnapshots())
           ys = sorted(getSpreaderCounts())
           xs = [1,2,3,5,7]
           ys = [1,2,3,5,10]
           plotfile = "{0}_{1}_{2}_{3}_{4}_{5}_{6}.png".format(infermode,curscore,smodel,tracefile,timefrac,noise,plottype)
           plotpath = "{0}/{1}".format(curplotfolder,plotfile)
           sideparams = [xs,ys,noise,timefrac]
         
        plotparams.append(plotpath)  
        if plottype in ["sprobparam","distparam","fracs","snapshots","spreadercount","noise"]:
           genPlot(algoparams,traceparams,plotparams,sideparams)
        elif plottype in ["timecount-frac","frac-timecount","frac-startcount","timecount-startcount"]:   
           genHeatMap(algoparams,traceparams,plotparams,sideparams)
        
if __name__ == "__main__":
    main()

