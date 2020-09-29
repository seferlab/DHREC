#include <stdio.h>
#include "graph.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <unordered_map>
using namespace std;
typedef Graph<int,int,int> GraphType;


GraphType readCutFile(char* filename){ 
    unordered_map <int, float> nodemap;
    unordered_map <tuple<int,int>, float> edgemap;
    
    #m["foo"] = 42;
    #cout << m["foo"] << endl;

    
    nodedict = {}
    string line;
    ifstream myfile (filename);
    if (myfile.is_open()) {
       while ( myfile.good() ) {
          getline (myfile,line);
          unsigned pos = str.find(" ");  
          int node = atoi(str.substr(0,pos).c_str());     
          float val = atof(str.substr(pos).c_str());
       }
       myfile.close();
    }
    else cout << "Unable to open file"; 
    
    GraphType *g = new GraphType(nodemap.size, edgemap.size);
    for(int curnode=0; curnode<nodemap.size; curnode++){
       g -> add_node(); 
       g -> add_tweights(curnode, /* capacities */ nodemap[curnode], nodemap[curnode]);          
    } 
    for (iter = edgemap.begin(); iter != edgemap.end(); ++iter) {
        g -> add_edge(iter[0], iter[1], /* capacities */  edgemap[iter], edgemap[iter] );  
    }

    return g
}
     

void WriteResult(g){
       
}



int main()  {
    char * filename = "";
    GraphType g = readCutFile(char* filename)
    int flow = g -> maxflow();
    
    printf("Flow = %d\n", flow);
    printf("Minimum cut:\n");
    if (g->what_segment(0) == GraphType::SOURCE)
       printf("node0 is in the SOURCE set\n");
    else
       printf("node0 is in the SINK set\n");
    if (g->what_segment(1) == GraphType::SOURCE)
       printf("node1 is in the SOURCE set\n");
    else
       printf("node1 is in the SINK set\n");

    delete g;

    return 0;
}
