import networkx as nx
import os
import sys
sys.path.append("./lib")
import myutilities as myutil


seens = set()
first = set()
second = set()
myfiles = set()
for mainfolder in myutil.listdirectories("."):
    if not mainfolder.startswith("result"):
       continue 
    #if mainfolder.find("si") == -1:
    #   continue 
    for file1 in myutil.listdirectories(mainfolder):
        path1 = "{0}/{1}".format(mainfolder,file1)
        for file2 in myutil.listdirectories(path1):
            path2 = "{0}/{1}".format(path1,file2)
            for file3 in myutil.listdirectories(path2):
                path3 = "{0}/{1}".format(path2,file3)
                for file4 in myutil.listdirectories(path3):
                    path4 = "{0}/{1}".format(path3,file4)
                    assert file4 in ["Spreader","History"]
                    #if file4 == "History":
                    #   continue
                    for file5 in myutil.listdirectories(path4):
                        path5 = "{0}/{1}".format(path4,file5)
                        seens.add(file5)
                        continue
                        #if file5.find("MatroidSub-False") != -1:
                        #   first.add(path5)
                        #   #code = "rm -rf {0}".format(path5)
                        #   #os.system(code)  
                        #elif file5.find("MatroidSub-continuous-False") != -1:
                        #   second.add(path5)     
                        if file5.find("RumorCentrality") != -1:
                           myfiles.add(path5)
                           code = "rm -rf {0}".format(path5)
                           os.system(code)  
                    
                        #if file5.find("Qsapmax") != -1:
                        #   code = "rm -rf {0}".format(path5)
                        #   os.system(code) 
                        #if file5.find("Qsap-reliability-False") != -1:
                        #   print file5
                        #   exit(1) 
                        #   newpath5 = path5.replace("Qsap","Qsapmin")
                        #   code = "mv {0} {1}".format(path5,newpath5)
                        #   os.system(code) 
                        #if file5.find("Qsapmin-reliability-False") != -1:
                        #   code = "rm -rf {0}".format(path5)
                        #   os.system(code)  
                        #   second.add(file5)   
                        
print seens                          
print len(first)
print len(second)
#print list(second)[0:4]
exit(1)                 
#set(['Qsapmax-False-greedy', 'NetSleuth', 'Independent-False', 'Pcdsvc-False', 'MatroidSub-search-False', 'Qsapmin-reliability-False', 'GreedySubCoverSingle-log-False-None', 'MinCut-False', 'Pcdsvc-True-5', 'MatroidSub-continuous-False', 'MatroidSub-False', 'GreedySubCoverSingle-log-True-None-3', 'GreedySubCoverSingle-log-True-None-2', 'GreedySubCoverSingle-log-True-None-5', 'RumorCentrality'])       
                                                                
