#include <stdio.h>
#include "graph.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <assert.h>
#include <map>
#include <set>

using namespace std;
typedef Graph<float,float,float> GraphType;

pair<map<int, char>,float> runFlow(const char* filename){ 
    map <int, float> s2node;
    map <int, float> node2t;
    map <pair<int,int>, pair<float,float> > edgemap;
    
    string line;
    bool edgeval = false;
    ifstream infile(filename);
    while(getline(infile, line)) {
        istringstream iss(line);
        if(iss.str() == "") {
          edgeval = true;
          continue;    
        }
        if(!edgeval){
          int node;float sval;float tval;
          if (!(iss >> node >> sval >> tval)){
	     break;
          }
          s2node[node] = sval; 
          node2t[node] = tval; 
        } 
        else{
          int node1;int node2;float value;float revvalue;
          if (!(iss >> node1 >> node2 >> value >> revvalue)){
	     break;
          }
          edgemap[make_pair(node1,node2)] = make_pair(value,revvalue);     
        }    
    }
    infile.close();
     
    map<int,int> node2index;
    map<int,float>::iterator nodeiter;   
    int index = 0; 
    for(nodeiter=s2node.begin(); nodeiter!=s2node.end(); ++nodeiter){ 
       int node = (*nodeiter).first;
       node2index[node] = index;
       index++;     
    }

    GraphType* g = new GraphType(s2node.size(), edgemap.size()); 
    for(int tnode=0; tnode < node2index.size() ; ++tnode){
        g -> add_node();
    }

    for(nodeiter=s2node.begin(); nodeiter!=s2node.end(); ++nodeiter){
       int node = (*nodeiter).first;  
       int nodeindex = node2index[node]; 
       g -> add_tweights(nodeindex,s2node[node],node2t[node]);
    }
     
    map< pair<int,int>, pair<float,float> >::iterator edgeiter;
    for (edgeiter = edgemap.begin(); edgeiter != edgemap.end(); edgeiter++) {
        int node1 = (*edgeiter).first.first;
	int node2 = (*edgeiter).first.second;
        int nodeindex1 = node2index[node1];
        int nodeindex2 = node2index[node2];
        g -> add_edge(nodeindex1,nodeindex2,edgemap[(*edgeiter).first].first,edgemap[(*edgeiter).first].second);  
    }
    
    map<int,int> index2node;
    map<int,int>::iterator tempiter;   
    int maxindex = -1000000;
    int minindex = 10000000;
    for(tempiter=node2index.begin(); tempiter!=node2index.end(); ++tempiter){ 
       int node = (*tempiter).first;
       int index = (*tempiter).second;
       index2node[index] = node;
       if(index > maxindex){
	  maxindex = index ;
       } 
       if(index < minindex){
	  minindex = index;
       } 
    }
    assert(minindex + maxindex == node2index.size()-1);   

    float flow = g -> maxflow();
    map<int, char> segment;
    int scount = 0, tcount = 0;
    for(int nodeindex = 0; nodeindex < maxindex+1; ++nodeindex){
       int node = index2node[nodeindex];
       if(g->what_segment(nodeindex) == GraphType::SOURCE){
	  scount += 1;  
	  segment[node] = 's';
       }
       else{
	  tcount += 1;
	  segment[node] = 't';
       }
    }
    printf("s-t result info\n");
    printf("%d\n",scount);
    printf("%d\n",tcount);
    assert(scount + tcount == segment.size());
    return make_pair(segment,flow);

}
    
 
void WriteRandomGraph(const char* filename){
        int nodecount = 10000;
        float p = 0.1;
        
        srand(time(NULL));
        
        ofstream myfile;
        myfile.open(filename);
        map<int,char> temp;
	for(int i=0; i<nodecount; i++){
            float w = rand() % 1000;
            w += ((double) rand() / (RAND_MAX));
            double r  = ((double) rand() / (RAND_MAX));
            if(r <= 0.5){ 
	       myfile << i << "\t" << 's' << "\t" << w << "\n"; 
               temp[i] = 's'; 
            }
            else{
               myfile << i << "\t" << 't' << "\t" << w << "\n"; 
               temp[i] = 't'; 
            }
        }
        myfile << "\n";
        for(int node1=0; node1<nodecount; node1++){
           for(int node2=node1+1; node2<nodecount; node2++){
               double r  = ((double) rand() / (RAND_MAX));
               if(r >= p){
		  continue;
               }
               float w = rand()%5;
               /*if(temp[node1] == temp[node2]){
                  w = 30000.0;
	       }*/
               w += ((double) rand() / (RAND_MAX));
               r  = ((double) rand() / (RAND_MAX));
               if(r <= 0.5){ 
	          myfile << node1 << "\t" << node2 << "\t" << w << "\t" << 0.0 << "\n";  
	       }
               else{
                  myfile << node1 << "\t" << node2 << "\t" << 0.0 << "\t" << w << "\n";  
               } 
           }
        }
        myfile.close();
  
}


void writeResult(map<int,char> outsegment, float flow, const char* outfilename){
     ofstream myfile;
     myfile.open(outfilename); 
     myfile << flow << "\n";
     map<int,char>::iterator segiter;
     for (segiter = outsegment.begin(); segiter != outsegment.end(); segiter++) {
         myfile << (*segiter).first << "\t" << outsegment[(*segiter).first] << "\n"; 
     } 
     myfile.close();    
}


int main(int argc, char* argv[]) {
    if(argc != 3 ){
       cout << "Not enough parameters\n";
       exit(0);
    }
    const char* filename = argv[1];
    const char* outfilename = argv[2];
    //WriteRandomGraph(filename);

    pair<map<int,char>,float> respair = runFlow(filename); 
    map<int,char> segment = respair.first;
    float flow = respair.second;
    printf("%d\n",segment.size());
    printf("%f\n",flow);
    writeResult(segment,flow,outfilename); 
    return 0;
   
}
