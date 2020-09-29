#include <stdio.h>
#include "graph.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
using namespace std;
typedef Graph<int,int,int> GraphType;

int main()
{
        //typedef Graph<int,int,int> GraphType;
	srand(time(NULL));

        int i,j;
        int nodecount = 30000;
	float p = 0.2;
        int edgecount = 220000;
        int ecount = 0;

        GraphType *g = new GraphType(/*estimated # of nodes*/ nodecount, /*estimated # of edges*/ edgecount); 

	for(i=0; i<nodecount; i++){
	    g -> add_node(); 
            int w = rand() % 50;
            double r  = ((double) rand() / (RAND_MAX));
            if(r <= 0.5){ 
               g -> add_tweights(i, /* capacities */  w, 0); 
            }
            else{ 
               g -> add_tweights(i, /* capacities */  0, w ); 
            }
        }

        int flag = 1;
        for(i=0; i<nodecount; i++){
	  for(j=i+1; j<nodecount; j++){ 
               int w = rand() % 50;
               double r  = ((double) rand() / (RAND_MAX));
               if(r <= 0.5){ 
	          g -> add_edge( i, j,    /* capacities */ w, 0 );
	       }
               else{
		 g -> add_edge( i, j,    /* capacities */ 0, w );
               } 
               ecount += 1; 
               if(ecount == edgecount){
                 flag = 0;
                 break; 
               } 
                   
           }
           if(flag == 0){
	     break;
	   }
	}
        
        printf("%d\n",ecount);
        printf("Starting flow algo\n");
	int flow = g -> maxflow();
	printf("Flow = %d\n", flow);
	printf("Minimum cut:\n");  
        map<pair<int,int>, char> segment;
        //map segment;
        int scount =0, tcount = 0;
        for(int node=0; node<nodecount ; node++){
          if(g->what_segment(node) == GraphType::SOURCE){
              scount += 1;  
              //segment[(node,node)] = 's';
              segment[make_pair(node,node)] = 's';
          }
          else{
	      tcount += 1;
              segment[make_pair(node,node)] = 't';
	      //segment[(node,node)] = 't';
          }
	  //cout << node << g->what_segment(node) << "\n";
        }
        printf("info\n");
	printf("%d\n",scount);
        printf("%d\n",tcount);
        printf("%d\n",segment.size());
        printf("%d\n",segment[make_pair(1,1)]);
        printf("%d\n",segment[make_pair(3,3)]);
        printf("%d\n",segment[make_pair(76,76)]);
        
        //my_map.find( key ) != my_map.end()
	return 0;
        
        
	/*if (g->what_segment(0) == GraphType::SOURCE)
	    printf("node0 is in the SOURCE set\n");
	else
	    printf("node0 is in the SINK set\n");
	if (g->what_segment(1) == GraphType::SOURCE)
	    printf("node1 is in the SOURCE set\n");
	else
	   printf("node1 is in the SINK set\n");*/

	delete g;

	return 0;
}
