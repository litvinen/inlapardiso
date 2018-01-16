// This file contains interface between sparse solver Pardiso and INLA library
// 14 January 2018
// KAUST, Bayesian Computational statistics and modeling, CEMSE, KAUST, 
// https://bayescomp.kaust.edu.sa/Pages/Home.aspx
// Dr. Alexander Litvinenko and Prof. Haavrd Rue

/* -------------------------------------------------------------------- */    
/* ..  We assume that the user is using the Fortran 1-based        */
/*     indexing and not 0-based C-notation.  */
/* -------------------------------------------------------------------- */ 

#include "graph-matrix-format.h"

/*The structure which contains all required data*/
struct data_storage{
  /* Internal solver memory pointer pt,                  */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
  /* or void *pt[64] should be OK on both architectures  */ 
  void    *pt[64]; 

  /* Pardiso control parameters. */
  int      iparm[64];
  double   dparm[64];
  int      maxfct, mnum, phase, error, msglvl, solver;

  /* Number of processors. */
  int      num_procs;
  int nnz;
  int n;
  int* perm;
  pspmatrix L;
};
typedef struct data_storage sdata_storage;
typedef sdata_storage* pdata_storage;
/*end*/

/*Sparse matrix structure */
struct spmatrix{
   int     nia = 0;
   int*    ia = NULL;
   int*    ja = NULL;
   double*  a = NULL;
   int      nnz = 0;
};
typedef struct spmatrix sspmatrix;
typedef sspmatrix* pspmatrix;
/*end*/


/* Sparse matrix structure as it is used in INLA*/
struct matrixQ{
  graph_t * g;
  Qfunc_t Q;
  void *arg;
};
typedef struct matrixQ smatrixQ;
typedef smatrixQ* pmatrixQ;
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
void alex_initialization(int flag)
{
   int      mtype = -2;        /* Real symmetric matrix */
      //       mtype = 1;   /* NOT real symmetric matrix*/      
   void    *pt[64];     /* Internal solver memory pointer pt,*/

      /* Pardiso control parameters. */
   int      iparm[64];
   double   dparm[64];

   pspmatrix Q = (pspmatrix) malloc(sizeof(sspmatrix));
   pspmatrix L = (pspmatrix) malloc(sizeof(sspmatrix));

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


//flag==0=="reordering"
//flag==1="symbolic factorization"
//flag==2="numerical factorization"
//flag==3="Cholesky factorization"
//flag==4="Partial inversion"
   if(flag==0) //0="reordering"
   {
      int      maxfct, mnum, phase, error, msglvl, solver;
      int      num_procs;/* Number of processors. */
      char    *var;        /* Auxiliary variables. */
      int      i, k;       /* Auxiliary variables. */
      double   ddum;              /* Double dummy */
      int      idum;              /* Integer dummy. */
     /* -------------------------------------------------------------------- */
     /* ..  Setup Pardiso control parameters.                                */
     /* -------------------------------------------------------------------- */
     error = 0;
     solver=0;/* use sparse direct solver */
    //iparm[12]=1;

   }
   if(flag==1) //1="symbolic factorization"
   {
   }
   if(flag==2) //2="numerical factorization"
   {}
   if(flag==3) //3="Cholesky factorization"
   {
      IPARM(33)=1;  //compute logdet

   }
   if(flag==4) //4="Partial inversion"
   {}

}  


void convert_F2C(pspmatrix Q)
{
/* -------------------------------------------------------------------- */    
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */ 
   int n=Q->n;
   int i=0;
    for (i = 0; i < n+1; i++) {
        Q->ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        Q->ja[i] -= 1;
    }
}


void convert_C2F(pspmatrix Q)
{
/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
   int n=Q->n;
   int i=0;
   for (i = 0; i < n+1; i++) 
   {
        Q->ia[i] += 1;
   }
   for (i = 0; i < nnz; i++) 
   {
        Q->ja[i] += 1;
   }
}

/* Read sparse matrix in matlab format */
void read_sparse_matrix( pspmatrix Q)
{
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
    Q->ia = n;
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
    Q->ia = ia;
    Q->ja = ja;
    Q->a  = a;
    Q->nnz=nnz;

    fclose( f5);
}   





// reordering, computed just once
int alex_reordering(pdata_storage mydata, pgraph_t mgraph, int* mypermutation)
{
    clock_t start, end;
    double cpu_time_used;

    alex_initialization(1); // initialization for reordering
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
    maxfct = 1;		/* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */
    pdata_storage mydata_storage = NULL;
/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
    phase = 11; //analysis
    read_sparse_matrix(mydata_storage);

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
int alex_symbolic_factorization(pdata_storage mystorage)
{
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
    phase = 11; //analysis
    read_sparse_matrix();
    a=mystorage->a;
    ia = mystorage->ia;
    ja = mystorage->ja;
    n = mystorage->n;
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

//Numerical Cholesky factorization
int alex_chol(pdata_storage mydata, pmatrix Q)
{
   int phase = 22;
   iparm[32] = 1; /* compute determinant */
    ini
   // start = clock();
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//  end = clock();
//  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
//  printf("!!!!!!!!!!!!Numerical factorization took %f seconds to execute \n", cpu_time_used ); 
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");
  return 0;
}




//Solve linear system Lx=b
int alex_solve_Lx(pdata_storage mydata,  double* b)
{

/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
    int phase = 33;

    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
   
//    start = clock();
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
//    end = clock();
//    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
//    printf("!!!!!!!!!!!!1Back substitution and iterative refinement took %f seconds to execute \n", cpu_time_used );
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }

    printf("\nSolve completed ... ");
    printf("\nThe solution of the system is: ");
//    for (i = 0; i < n; i++) {
//        printf("\n x [%d] = % f", i, x[i] );
//    }
    printf ("\n\n");


  return 0;
}


//Solve linear system L^Tx=b
int alex_solve_LTx(pdata_storage mydata,  double* b)
{
  return 0;
}


//Solve linear system LL^Tx=b
int alex_solve_LLTx(pdata_storage mydata,  double* b)
{
  alex_solve_Lx(mydata,  b);
  alex_solve_LTx(mydata, b);
    
  return 0;
}

//Computing partial inverse, means solving linear system foer specific RHS
int alex_inv(pdata_storage mydata,  pmatrix Q)
{
  /*On entry: IPARM(36) will control the selected inversion process based on the internal L and
  U factors. If IPARM(36) = 0 and PHASE = -22, PARDISO will overwrite these factors with
  the selected inverse elements of A −1 . If IPARM(36) = 1 and PHASE = -22, PARDISO will
  not overwrite these factors. It will instead allocate additional memory of size of the numbers of
  elements in L and U to store these inverse elements.*/
  return 0;
}

//Compute log-det
double alex_log_det(data_storage* mydata)
{
   double logdet=0.0;
   IPARM(33)=1;  //compute logdet
/*
  IPARM (33) — Determinant of a matrix.
  Input
  On entry: IPARM(33)=1 will compute the determinant of matrices and will return in DPARM(33)
  the real part of the determinant and, if necessary, in DPARM(32) the complex part of the deter-
  minant.
  On output: The parameter returns the natural logarithm of the determinant of a sparse matrix
  A. The default value of IPARM(33) is 0.
*/
  return 0;
}


//Finalize and deallocate memory

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */    
int alex_finalize(pdata_storage mydata)
{
    int phase = -1;                 /* Release internal memory. */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

}

int alex_clean_mydata(pdata_storage  mydata)
{
   /* To clean what is necessary in mydata->... */
}


/* To test the code run this procedure */
int alex_main()
{
   clock_t start, end;
   double cpu_time_used;
   pdata_storage mydata;

   int nnz;
   int n;
   int* perm;
   pspmatrix L;
   int mtype = -2;        /* Real symmetric matrix */

    
   double* b; /* RHS and solution vectors. */
   double* x;
   int nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    void    *pt[64]; 
    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;
    int      num_procs;/* Number of processors. */
    /* Auxiliary variables. */
    char    *var;
    int      i, k, j;
    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */
    pspmatrix Q = NULL, Qinv=NULL, L=NULL;
    pgraph_t mgraph = NULL;
    int* mypermutation = NULL;

    double* b = NULL;     /* RHS and solution vectors. */
    double* x = NULL;
    int*  mypermutation = NULL;

    mydata->mtype = -2;
    mydata->rhs = 1;

    mypermutation = (int*)malloc(n*sizeof(int)) ;
    for( j = 0; j < n; j++)
      mypermutation[j] = j;







    read_sparse_matrix(Q, mydata);
    convert_C2F(Q);
    alex_initialization(0, mydata); //flag==0=="reordering"
    alex_reordering(mydata, mgraph, mypermutation);
    alex_clean_mydata(mydata);

    alex_initialization(1, mydata); //flag==1="symbolic factorization"
    alex_initialization(2, mydata); //flag==2="numerical factorization"
    alex_symbolic_factorization(mydata);
    alex_clean_mydata(mydata);

    alex_initialization(3, mydata); //flag==3="Cholesky factorization"
    alex_chol(mydata, Q , L);
    alex_clean_mydata(mydata);

    alex_initialization(5, mydata); //flag==5="Solution"
    alex_solve_LLTx(mydata, b, x);
    alex_clean_mydata(mydata);

    alex_initialization(4, mydata); //flag==4="Partial inversion"
    alex_inv(mydata, Q, Qinv);
    alex_clean_mydata(mydata);
    convert_F2C(Q);

}


