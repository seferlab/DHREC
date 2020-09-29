#Cut utilities
import networkx as nx
import sys
import os


def storeBoykovGraph(cutG,s,t,filename):
    """stores graph in format for boykov algo
    Args:
       cutG:
       s:
       t:
       filename:
    Returns:
    """
    with open(filename,"w") as outfile:
        for node in cutG.nodes():
            if node not in [s,t]:
               sval,tval = 0.0, 0.0
               if cutG.has_edge(s,node):
                  sval = cutG[s][node]["weight"] 
               if cutG.has_edge(node,t):
                  tval = cutG[node][t]["weight"]    
               outfile.write("{0}\t{1}\t{2}\n".format(node,sval,tval))
        outfile.write("\n")        
        seenedges = set()    
        for node1,node2 in cutG.edges():
            prenode, afternode = node1,node2
            if node2 < node1:
               prenode, afternode = node2,node1
            if (prenode,afternode) in seenedges:
               continue
            seenedges.add((prenode,afternode))
            forflow,revflow = 0.0,0.0
            if prenode not in [s,t] and afternode not in [s,t]:
               if cutG.has_edge(prenode,afternode):
                  forflow = cutG[prenode][afternode]["weight"] 
               if cutG.has_edge(afternode,prenode):
                  revflow = cutG[afternode][prenode]["weight"] 
               outfile.write("{0}\t{1}\t{2}\t{3}\n".format(prenode,afternode,forflow,revflow))
