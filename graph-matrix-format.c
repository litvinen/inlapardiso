// Prof. Haavard Rue
// A framework to work with sparse matrices.
// haavard.rue@kaust.edu.sa
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "graph-matrix-format.h"
//#include <boost/concept_check.hpp>



// this example shows the internal graph structure and the function defining the matrix.



// the full binary matrix is stored here
/*typedef struct {
	int n;						       // size
	int *nnbs;					       // number of neigbours
	int **nbs;					       // list of neighbours
} graph_t;
*/
// this functions returns element Q_ij for i=j or i~j
/*typedef double Qfunc_t(int i, int j, void *arg);
graph_t *read_graph(char *filename);
int print_graph(graph_t * g);
int print_Q(graph_t * g, Qfunc_t Q, void *arg);
double Qdemo(int i, int j, void *arg);
*/

graph_t *read_graph(char *filename)
{
#define READ_INT(location) fscanf(fp, "%d", &(location))

	int i, k, kk;
	graph_t *g = (graph_t *) calloc(1, sizeof(graph_t));
	FILE *fp = fopen(filename, "r");

	READ_INT(g->n);
	g->nnbs = (int *) calloc(g->n, sizeof(int));
	g->nbs = (int **) calloc(g->n, sizeof(int *));

	for (k = 0; k < g->n; k++) {
		READ_INT(i);
		READ_INT(g->nnbs[i]);
		g->nbs[i] = (int *) calloc(g->nnbs[i], sizeof(int));
		for (kk = 0; kk < g->nnbs[i]; kk++) {
			READ_INT(g->nbs[i][kk]);
		}
	}
#undef READ_INT
	return (g);
}

int print_graph(graph_t * g)
{
	int i, j;
	printf("n = %d\n", g->n);
	for (i = 0; i < g->n; i++) {
		printf("node %d has %d neighbours\t", i, g->nnbs[i]);
		for (j = 0; j < g->nnbs[i]; j++)
			printf(" %d", g->nbs[i][j]);
		printf("\n");
	}
	return (0);
}

int print_Q(graph_t * g, Qfunc_t Q, void *arg)
{
	int i, j, jj;
	printf("dim(Q) = %d\n", g->n);
	for (i = 0; i < g->n; i++) {
		// note that the diagonal are not neibours...
		printf("Q(%1d, %1d) = %g\n", i, i, Q(i, i, arg));
		for (jj = 0; jj < g->nnbs[i]; jj++) {
            j = g->nbs[i][jj];
			printf("\tQ(%1d, %1d) = %g\n", i, j, Q(i, j, arg));
		}
	}
	return (0);
}

double Qdemo(int i, int j, void *arg)
{
    // an example of a ``matrix''
    graph_t *g = (graph_t *) arg;
    return (i == j ? g->nnbs[i] + 1 : -1);
}




void convert2CSR(psspmatrix S, graph_t * g, Qfunc_t Q, void *arg)
{
  int n = 0, nnz=0;
  int i, jj, j;
  int k = 0;
  
  printf("rows = %d, cols = %d\n", g->n, g->n);
  for (i = 0; i < g->n; i++) 
     nnz=nnz + g->nnbs[i];
  
  S->nnz = nnz;
  S->n = g->n; 
  S->rows = g->n; 
  S->a = (double*)malloc(nnz * sizeof(double));
  S->ja =(int*)malloc(nnz*sizeof(int));
  S->ia=(int*)malloc(g->n*sizeof(int));
  S->ia[0]=0;
  
  for( k = 1; k <= g->n; k++)
  {
      S->ia[k] = S->ia[k-1] + g->nnbs[k-1]; 
  }
  for (i = 0; i < g->n; i++) {
    for (jj = 0; jj < g->nnbs[i]; jj++)
    {    
        j = g->nbs[i][jj]-1;
        S->a[k] = Q(i, j, arg);
        k++;
    }
  }

  k=0;
  for (i = 0; i < g->n; i++) {
    for (j = 0; j < g->nnbs[i]; j++)
    {    
        S->ja[k] = g->nbs[i][j]-1;
        k++;
    }
  }
  
}


void convert2CSR_alex(psspmatrix S, graph_t * g, double* values)
{
  int n = 0, nnz=0;
  int i, jj, j;
  int k = 0;
  
  printf("rows = %d, cols = %d\n", g->n, g->n);
  for (i = 0; i < g->n; i++) 
     nnz=nnz + g->nnbs[i];
  
  S->nnz = nnz;
  S->n = g->n; 
  S->rows = g->n; 
  S->a = (double*)malloc(nnz * sizeof(double));
  S->ja =(int*)malloc(nnz*sizeof(int));
  S->ia=(int*)malloc(g->n*sizeof(int));
  S->ia[0]=0;
  
  for( k = 1; k <= g->n; k++)
     S->ia[k] = S->ia[k-1] + g->nnbs[k-1]; 
  for (i = 0; i < S->nnz; i++) 
     S->a[i] = values[i];
  
  k=0;
  for (i = 0; i < g->n; i++) {
    for (j = 0; j < g->nnbs[i]; j++)
    {    
        S->ja[k] = g->nbs[i][j]-1;
        k++;
    }
  }
  
}
void print_CSR(psspmatrix S)
{
    int k=0;
    printf("Array A: \n"); 
    for( k = 0; k < S->nnz; k++)
       printf("%3.3g, ", S->a[k]); 
    printf(" \n"); 
    for( k = 0; k <= S->rows; k++)
       printf("%d, ", S->ia[k]); 
    printf(" \n"); 
    for( k = 0; k < S->nnz; k++)
       printf("%d, ", S->ja[k]); 
    printf(" \n"); 
  
}

int main(int argc, char **argv)
{
	graph_t *g;
    psspmatrix S;
    double* v;
    
    v=(double* )malloc(4*sizeof(double));
    v[0]=5;
    v[1]=8;
    v[2]=3;
    v[3]=6;
    
    S=(psspmatrix )malloc(sizeof(sspmatrix ));
    g = read_graph("germany.graph.txt");
   // g = read_graph("minitest.graph.txt");
    /* CRS matrix is
     0 0 0 0
     5 8 0 0
     0 0 3 0
     0 6 0 0
     */
    print_graph(g);
	print_Q(g, Qdemo, (void *) g);
    convert2CSR(S, g, Qdemo, (void *) g);
//    convert2CSR_alex(S, g, v);
    print_CSR(S);
    /*
     should obtain
     A= [5 8 3 6]
     IA=[0 0 2 3 4]
     JA=[0 1 2 1]
     */
    
	return (0);
}
