// Prof. Haavard Rue
// A framework to work with sparse matrices.
// haavard.rue@kaust.edu.sa
// this example shows the internal graph structure and the function defining the matrix.
#ifndef _MYGRAPH
#define _MYGRAPH

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

// the full binary matrix is stored here
typedef struct {
	int n;						       // size
	int *nnbs;					       // number of neigbours
	int **nbs;					       // list of neighbours
} graph_t;


typedef struct graph_t sgraph_t;
typedef sgraph_t* psgraph_t;


// this functions returns element Q_ij for i=j or i~j
typedef double Qfunc_t(int i, int j, void *arg);
graph_t *read_graph(char *filename);
int print_graph(graph_t * g);
//int print_Q(graph_t * g, Qfunc_t Q, void *arg);
int print_Q(psgraph_t  g,  void *arg);
double Qdemo(int i, int j, void *arg);


#endif
