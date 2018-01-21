// This file contains interface functions between the sparse solver Pardiso and the INLA library
// 14 January 2018
// KAUST, Bayesian Computational statistics and modeling, CEMSE, KAUST, 
// https://bayescomp.kaust.edu.sa/Pages/Home.aspx
// Dr. Alexander Litvinenko and Prof. Haavard Rue
#ifndef _MYINTERFACE
#define _MYINTERFACE

/* -------------------------------------------------------------------- */    
/* ..  We assume that the user is using the Fortran 1-based        */
/*     indexing and not 0-based C-notation.  */
/* -------------------------------------------------------------------- */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#include "graph-matrix-format.h"

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

/*The structure which contains all required data*/
struct data_storage{
  /* Internal solver memory pointer pt,                  */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
  /* or void *pt[64] should be OK on both architectures  */ 
  void    *pt[64]; 

  /* Pardiso control parameters. */
  int      iparm[64];
  double   dparm[64];
  int  mtype;
  int  nrhs;

  int      maxfct;
  int      mnum;
  int      phase;
  int      error;
  int      msglvl;
  int      solver;
  int      num_procs;   /* Number of processors. */
  int nnz;
  int n;
  int* perm;
  char    *var;        /* Auxiliary variables. */
  int      i, k;       /* Auxiliary variables. */
  double   ddum;              /* Double dummy */
  int      idum;              /* Integer dummy. */
  //pspmatrix Q = NULL;
  //pspmatrix L = NULL;
  psspmatrix L;
  psspmatrix Q;
  psspmatrix Qinv;
  //psspmatrix Q = (psspmatrix) malloc(sizeof(sspmatrix));

};
typedef struct data_storage sdata_storage;
typedef sdata_storage* pdata_storage;
/*end*/


/* Sparse matrix structure as it is used in INLA*/
struct INLA_mtx{
  psgraph_t  g;
  //Qfunc_t Q;
  void *arg;
};
typedef struct INLA_mtx sINLA_mtx;
typedef sINLA_mtx* pINLA_mtx;
/*end*/



//struct spvector{
//  int n;
//  int* ind;
//  double* vals;
//};
//typedef struct spvector sspvector;
//typedef sspvector* pspvector;



//typedef struct {
//	int n;						       // size
//	int *nnbs;					       // number of neigbours
//	int **nbs;					       // list of neighbours
//} graph_t;
//typedef struct graph_t sgraph_t;
//typedef sgraph_t* pgraph_t;


/*Coincidence matrix (GRAPH)*/
//typedef struct coin_matrix scoin_matrix;
//typedef scoin_matrix* pcoin_matrix;

/*The structure which contains all required data*/
//struct coin_matrix{
//  int rows;
//  int cols;
//  int* ia;
//  int* ja;
//};
//int*    ia = NULL;
//int*    ja = NULL;
//double*  a = NULL;


/* Allocation and initialization of PARDISO variables*/
int alex_initialization(int flag, pdata_storage mydata);


void convert_F2C(psspmatrix Q);


void convert_C2F(psspmatrix Q);

/* Read sparse matrix in CSR format */
void read_sparse_matrix( psspmatrix Q);

// reordering, computed just once
int alex_reordering(pdata_storage mydata, psgraph_t mgraph, int* mypermutation);

//compute the symbolic factorization, computed just once
int alex_symbolic_factorization(pdata_storage mydata);

//Numerical Cholesky factorization
//int alex_chol(pdata_storage mydata, pINLA_mtx Q);
int alex_chol(pdata_storage mydata, psspmatrix Q, psspmatrix L);

//Solve linear system Lx=b
int alex_solve_Lx(pdata_storage mydata,  double* b,  double* x);


//Solve linear system L^Tx=b
int alex_solve_LTx(pdata_storage mydata,  double* b,  double* x);


//Solve linear system LL^Tx=b
int alex_solve_LLTx(pdata_storage mydata,  double* b,  double* x);

//Computing partial inverse, means solving linear system foer specific RHS
//int alex_inv(pdata_storage mydata,  pINLA_mtx Q,  pINLA_mtx Qinv);
int alex_inv(pdata_storage mydata,  psspmatrix Q,  psspmatrix Qinv);

//Compute log-det
double alex_log_det(pdata_storage mydata);


//Finalize and deallocate memory

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */    
int alex_finalize(pdata_storage mydata);

int alex_clean_mydata(pdata_storage  mydata);

/* To test the code run this procedure */
void init_mydata(pdata_storage mydata);

int main();
#endif
