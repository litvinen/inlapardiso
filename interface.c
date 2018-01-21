// This file contains interface functions between the sparse solver Pardiso and the INLA library
// 14 January 2018
// KAUST, Bayesian Computational statistics and modeling, CEMSE, KAUST, 
// https://bayescomp.kaust.edu.sa/Pages/Home.aspx
// Dr. Alexander Litvinenko and Prof. Haavard Rue

/* -------------------------------------------------------------------- */    
/* ..  We assume that the user is using the Fortran 1-based        */
/*     indexing and not 0-based C-notation.  */
/* -------------------------------------------------------------------- */ 

/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
*/
#include "interface.h"
//#include <boost/iterator/iterator_concepts.hpp>
//#include <boost/concept_check.hpp>


/* Allocation and initialization of PARDISO variables*/
int alex_initialization(int flag, pdata_storage mydata)
{
  int num_procs=0;
  char    *var;
    
//flag==-1=="initial/common initialization of PARDISO"
//flag==0=="reordering"
//flag==1="symbolic factorization"
//flag==2="numerical factorization"
//flag==3="Cholesky factorization"
//flag==4="Partial inversion"


   if(flag==-1) //general initialization
   { 
     pardisoinit (mydata->pt,  &(mydata->mtype), &(mydata->solver), mydata->iparm, mydata->dparm, 
        &(mydata->error)); 
     if (mydata->error != 0) 
     {
        if (mydata->error == -10 )
           printf("No license file found \n");
        if (mydata->error == -11 )
           printf("License is expired \n");
        if (mydata->error == -12 )
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
     mydata->iparm[2]  = num_procs;

     mydata->maxfct = 1;		/* Maximum number of numerical factorizations.  */
     mydata->mnum   = 1;         /* Which factorization to use. */
     mydata->msglvl = 1;         /* Print statistical information  */
     mydata->error  = 0;         /* Initialize error flag */



     pardiso_chkmatrix  (&(mydata->mtype), &(mydata->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->error));
     if (mydata->error != 0) {
        printf("\nERROR in consistency of matrix: %d", mydata->error);
        exit(1);
     }
   }


   if(flag==0) //0="reordering"
   {
     /* -------------------------------------------------------------------- */
     /* ..  Setup Pardiso control parameters.                                */
     /* -------------------------------------------------------------------- */
      mydata->phase = 11; 
      mydata->error = 0;
    //iparm[12]=1;

   }
   if(flag==1) //1="symbolic factorization"
   {
      mydata->phase = 11;
   }
   if(flag==2) //2="numerical factorization"
   {
      mydata->phase = 22;
      mydata->iparm[32] = 1; /* compute determinant */
   }
   if(flag==3) //3="Cholesky factorization"
   {
      mydata->phase = 12;
   }
   if(flag==4) //4="Selected inversion"
   {
      mydata->phase = -22;
      mydata->iparm[35]  = 1; /*  no not overwrite internal factor L */ 

   }
   if(flag==5) //4="Solve for x"
   {
      mydata->iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
      mydata->phase = 33;
   }

   return 1;
}  


void convert_F2C(psspmatrix Q)
{
/* -------------------------------------------------------------------- */    
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */ 
   int n=Q->n;
   int nnz=Q->nnz;
   int i=0;

   for (i = 0; i < n+1; i++) {
        Q->ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        Q->ja[i] -= 1;
    }
}


void convert_C2F(psspmatrix Q)
{
/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
   int n=Q->n;
   int nnz=Q->nnz;
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
void read_sparse_matrix( psspmatrix Q)
{
    int    n = 0;
    int    nnz = 0;
    FILE *f5;
    int anzv=0, j=0, row_index=0, col_index=0;
    double elem=0.0;
    double* a=NULL;
    int* ia=NULL;
    int* ja=NULL;
    
    f5 = fopen( "ia.txt", "r");
    if( f5 == NULL)
    {
      printf( "Can't open the file\n");
      exit (1);
    }
    fscanf( f5, "%d\n", &anzv);
    n = anzv-1;
    Q->n=n;
    Q->ia=(int*)malloc((n+1)*sizeof(int));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%d", &row_index);
      ia[j] = row_index; 
    }
    fclose( f5);
    nnz = ia[n];
    Q->nnz=nnz;

    f5 = fopen( "ja.txt", "r");
    if( f5 == NULL)
    {
      printf( "Can't open the file\n");
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
      printf( "Can't open the file\n");
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
    
    fclose( f5);
}   





// reordering, computed just once
int alex_reordering(pdata_storage mydata, psgraph_t mgraph, int* mypermutation)
{
    clock_t start, end;
    double cpu_time_used;

    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
	     &(mydata->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error), mydata->dparm);
    
    if (mydata->error != 0) {
        printf("\nERROR during symbolic factorization: %d", mydata->error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", mydata->iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", mydata->iparm[18]);

  return 0;
}




//compute the symbolic factorization, computed just once
int alex_symbolic_factorization(pdata_storage mydata)
{
    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
	     &(mydata->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error), mydata->dparm);
    
    if (mydata->error != 0) {
        printf("\nERROR during symbolic factorization: %d", mydata->error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", mydata->iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", mydata->iparm[18]);

  return 0;
}

//Numerical Cholesky factorization
//int alex_chol(pdata_storage mydata, pINLA_mtx Q)
int alex_chol(pdata_storage mydata, psspmatrix Q, psspmatrix L)
{
    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error),  mydata->dparm);
    if (mydata->error != 0) {
        printf("\nERROR during numerical factorization: %d", mydata->error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");
  return 0;
}




//Solve linear system Lx=b
int alex_solve_Lx(pdata_storage mydata,  double* b,  double* x)
{

/* ..  Back substitution and iterative refinement.                      */
    pardiso ( mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), b, x, &(mydata->error),  mydata->dparm);
    if (mydata->error != 0) {
        printf("\nERROR during solution: %d", mydata->error);
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
int alex_solve_LTx(pdata_storage mydata,  double* b,  double* x)
{
    clock_t start, end;
    double cpu_time_used;
    int i=0;
    int n = mydata->n;
/* -------------------------------------------------------------------- */    
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
    start = clock();
    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), b, x, &(mydata->error),  mydata->dparm);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("!!!!!!!!!!!!1Back substitution and iterative refinement took %f seconds to execute \n", cpu_time_used );
    if (mydata->error != 0) {
        printf("\nERROR during solution: %d", mydata->error);
        exit(3);
    }

    printf("\nSolve completed ... ");
    printf("\nThe solution of the system is: ");
    for (i = 0; i < n; i++) {
//        printf("\n x [%d] = % f", i, x[i] );
    }
    printf ("\n\n");

  return 0;
}


//Solve linear system LL^Tx=b
int alex_solve_LLTx(pdata_storage mydata,  double* b,  double* x)
{
  alex_solve_Lx(mydata,  b, x);
  alex_solve_LTx(mydata, b, x);
    
  return 0;
}

//Computing partial inverse, means solving linear system foer specific RHS
int alex_inv(pdata_storage mydata,   psspmatrix Q,   psspmatrix Qinv)
{
  /*On entry: IPARM(36) will control the selected inversion process based on the internal L and
  U factors. If IPARM(36) = 0 and PHASE = -22, PARDISO will overwrite these factors with
  the selected inverse elements of A −1 . If IPARM(36) = 1 and PHASE = -22, PARDISO will
  not overwrite these factors. It will instead allocate additional memory of size of the numbers of
  elements in L and U to store these inverse elements.*/
  int k=0, i=0;
  double* x;
  double* b;
  for (i = 0; i < Q->n; i++) {
      b[i] = i;
      x[i] = 0.0;
  }

   if (mydata->solver == 0)
   {
    	printf("\nCompute Diagonal Elements of the inverse of A ... \n");
	    mydata->phase = -22;
        mydata->iparm[35]  = 1; /*  no not overwrite internal factor L */ 

  
        pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase), &(mydata->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), b, x, &(mydata->error),  mydata->dparm);

//        printf("!!!!!!!!!!!!1Inverse factorization took %f seconds to execute \n", cpu_time_used );
 /* print diagonal elements */
       for (k = 0; k < mydata->n; k++)
       {
            int j = mydata->Q->ia[k]-1;      
            printf ("Diagonal element of A^{-1} = %d %d %32.24e\n", k, mydata->Q->ja[j]-1, mydata->Q->a[j]);
       }

   } 


  return 0;
}

//Compute log-det
double alex_log_det(pdata_storage mydata)
{
   double logdet=0.0;

   pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error),  mydata->dparm);
/*
  IPARM (33) — Determinant of a matrix.
  Input
  On entry: IPARM(33)=1 will compute the determinant of matrices and will return in DPARM(33)
  the real part of the determinant and, if necessary, in DPARM(32) the complex part of the deter-
  minant.
  On output: The parameter returns the natural logarithm of the determinant of a sparse matrix
  A. The default value of IPARM(33) is 0.
*/

  mydata->Q->logdet=   mydata->dparm[33];
  return mydata->dparm[33];
}


//Finalize and deallocate memory

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */    
int alex_finalize(pdata_storage mydata)
{
    int phase = -1;                 /* Release internal memory. */
    
    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->n), &(mydata->ddum), mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error),  mydata->dparm);

}

int alex_clean_mydata(pdata_storage  mydata)
{
   /* To clean what is necessary in mydata->... */
}


void init_mydata(pdata_storage mydata)
{
    
   mydata->Q = (psspmatrix) malloc(sizeof(sspmatrix));
   mydata->L = (psspmatrix) malloc(sizeof(sspmatrix));
   mydata->Qinv = (psspmatrix) malloc(sizeof(sspmatrix));
    
}


/* To test the code run this procedure */
int main()
{
   clock_t start, end;
   double cpu_time_used;
   pdata_storage mydata= NULL;
  

   int nnz;
   int n;
   int* perm;
   
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
    psspmatrix Q = NULL, Qinv=NULL, L=NULL;
    psgraph_t mgraph = NULL;
    int* mypermutation = NULL;

    solver=0;/* use sparse direct solver */
    mtype=-2;
    nrhs=1;
    mydata= (pdata_storage)malloc(sizeof(sdata_storage ));
  ;
printf("flag -11\n");
    init_mydata(mydata);
printf("flag 0\n");
    
    mydata->mtype = mtype;
    mydata->nrhs = nrhs;
    mydata->solver = solver;
printf("flag 1\n");
    read_sparse_matrix(Q);
printf("flag 2\n");
    n=Q->n;
    mypermutation = (int*)malloc(n*sizeof(int)) ;
    x = (double*)malloc(n*sizeof(double)) ;
    b = (double*)malloc(n*sizeof(double)) ;
    for( j = 0; j < n; j++)
    {
      mypermutation[j] = j;
      x[i]=0.0; 
      b[i]=1.0; 
    } 



    convert_C2F(Q);
printf("flag 3\n");


    alex_initialization(0, mydata); //flag==0=="reordering"
    alex_reordering(mydata, mgraph, mypermutation);
    alex_clean_mydata(mydata);

    alex_initialization(1, mydata); //flag==1="symbolic factorization"
    alex_initialization(2, mydata); //flag==2="numerical factorization"
    alex_symbolic_factorization(mydata);
    alex_clean_mydata(mydata);

    alex_initialization(3, mydata); //flag==3="Cholesky factorization"
    alex_chol(mydata, Q , L);  //also computes log determinant
    alex_clean_mydata(mydata);

    alex_initialization(5, mydata); //flag==5="Solution"
    alex_solve_LLTx(mydata, b, x);
    alex_clean_mydata(mydata);

    alex_initialization(4, mydata); //flag==4="Partial inversion"
    alex_inv(mydata, Q, Qinv);
    alex_clean_mydata(mydata);
    convert_F2C(Q);
    alex_finalize(mydata);
    if(x!=NULL) free(x);
    if(b!=NULL) free(b);
    
    
    return 0;
}


