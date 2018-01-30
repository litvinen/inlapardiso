// Prof. Haavard Rue
// A framework to work with sparse matrices.
// haavard.rue@kaust.edu.sa
// this example shows the internal graph structure and the function defining the matrix.
#ifndef _MYGRAPH
#define _MYGRAPH

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
/*Sparse matrix structure */


/*Sparse matrix structure */
struct spmatrix{
   int     n ;
   int     nia ;
   int     cols;
   int     rows;
   int*    ia;
   int*    ja;
   double*  a ;
   int      nnz;
   int      logdet;
};
typedef struct spmatrix sspmatrix;
typedef sspmatrix* psspmatrix;
/*end*/

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
double Qdemo(int i, int j, void *arg);
int print_Q(graph_t * g, Qfunc_t Q, void *arg);
//void convert2CSR(graph_t * g, psspmatrix S);
void convert2CSR(psspmatrix S, graph_t * g, Qfunc_t Q, void *arg);
void convert2CSR_alex(psspmatrix S, graph_t * g, double* values);
void print_CSR(psspmatrix S);

#endif
