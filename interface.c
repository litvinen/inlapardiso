// This file contains interface between sparse solver Pardiso and INLA library
// 14 January 2018
// KAUST, Bayesian Computational statistics and modeling, CEMSE, KAUST, 
// https://bayescomp.kaust.edu.sa/Pages/Home.aspx
// Dr. Alexander Litvinenko and Prof. Haavrd Rue

#include "graph-matrix-format.h"

/*The structure which contains all required data*/
struct data_storage{
  int rows;
  int cols;
  int nnz;
  int* perm;
  int* ia;
  int* ja;
  double* vals;
};
typedef struct data_storage sdata_storage;
typedef sdata_storage* pdata_storage;



struct spvector{
  int n;
  int* ind;
  double* vals;
};
typedef struct spvector sspvector;
typedef sspvector* pspvector;

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
int*    ia = NULL;
int*    ja = NULL;
double*  a = NULL;
    int      nnz = 0;


void read_sparse_matrix()
{
    /* Matrix data. */
    int    n = 0;
    FILE *f5;
    int anzv=0, j=0, row_index=0, col_index=0;
    
    f5 = fopen( "ia.txt", "r");
    if( f5 == NULL)
    {
      printf( "Datei konnte nicht geoeffnet werden\n");
      exit (1);
    }
    fscanf( f5, "%d\n", &anzv);
    n = anzv-1;
    ia=(int*)malloc((n+1)*sizeof(int));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%d", &row_index);
      ia[j] = row_index; 
    }
    fclose( f5);
    nnz = ia[n];
    f5 = fopen( "ja.txt", "r");
    if( f5 == NULL)
    {
      printf( "Datei konnte nicht geoeffnet werden\n");
      exit (1);
    }
    fscanf( f5, "%d\n", &anzv); 
    ja=(int*)malloc((anzv+1)*sizeof(int));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%d", &col_index);
      ja[j] = col_index; 
    }
    fclose( f5);
    f5 = fopen( "a.txt", "r");
    if( f5 == NULL)
    {
      printf( "Datei konnte nicht geoeffnet werden\n");
      exit (1);
    }
    fscanf( f5, "%d\n", &anzv); 
    a=(double*)malloc((anzv+1)*sizeof(double));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%lf", &elem);
      a[j] = elem; 
    }
    fclose( f5);
}   


// reordering, computed just once
int alex_reordering(pdata_storage mydata_storage, pgraph_t mgraph)
{
    clock_t start, end;
    double cpu_time_used;


    int      mtype = -2;        /* Real symmetric matrix */
    //    int      mtype = 1;   /* NOT real symmetric matrix*/      
    void    *pt[64];     /* Internal solver memory pointer pt,*/

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;
    int      num_procs;/* Number of processors. */
    char    *var;        /* Auxiliary variables. */
    int      i, k;       /* Auxiliary variables. */
    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */
    FILE *f5;
    int anzv=0, j=0, row_index=0, col_index=0;
/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters.                                */
/* -------------------------------------------------------------------- */
    error = 0;
    solver=0;/* use sparse direct solver */
    //iparm[12]=1;
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error); 
    if (error != 0) 
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         return 1; 
    }
    else
        printf("[PARDISO]: License check was successful ... \n");
    
    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;
    maxfct = 1;		/* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
    phase = 11; 
    read_sparse_matrix();

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
    
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
  return 0;
}

//compute the symbolic factorization, computed just once
int alex_symbolic_factorisation(pdata_storage mydata_storage)
{
  return 0;
}

//Cholesky factorization
int alex_chol(pdata_storage mydata_storage, pmatrix Q)
{
  return 0;
}

//Solve linear syatem Lx=b
int alex_solve_Lx(pdata_storage mydata_storage,  double* b)
{
  return 0;
}


//Solve linear system L^Tx=b
int alex_solve_LTx(pdata_storage mydata_storage,  double* b)
{
  return 0;
}


//Solve linear system LL^Tx=b
int alex_solve_LLTx(pdata_storage mydata_storage,  double* b)
{
  return 0;
}

//Computing partial inverse, means solving linear system foer specific RHS
int alex_inv(pdata_storage mydata_storage,  pmatrix Q)
{
  return 0;
}

//Compute log-det
int alex_log_det(data_storage* mydata_storage)
{
  return 0;
}





