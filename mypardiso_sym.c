/* -------------------------------------------------------------------- */
/*      Example program to show the use of the "PARDISO" routine        */
/*      on symmetric linear systems                                     */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Olaf Schenk, Institute of Computational Science                 */
/*      Universita della Svizzera italiana, Lugano, Switzerland.        */
/*      Email: olaf.schenk@usi.ch                                       */
/* -------------------------------------------------------------------- */
//kw13919:litvina~/PARDISO> gcc mypardiso_unsym.c -g -o mypardiso_unsym.exe -L/home/litvina/PARDISO libpardiso500-GNU461-X86-64.so -llapack -lblas  -lgfortran -fopenmp -lpthread -lm


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/* PARDISO prototype. */
void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);


int main( void ) 
{
    /* Matrix data. */
    int    n = 0;
    int*    ia = NULL;
    int*    ja = NULL;
    double*  a = NULL;
    double elem=0.0;

    clock_t start, end;
    double cpu_time_used;


    int      nnz = 0;
    int      mtype = -2;        /* Real symmetric matrix */
//    int      mtype = 1;       

    /* RHS and solution vectors. */
    double* b;
    double* x;
    int      nrhs = 1;          /* Number of right hand sides. */

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

    /* Auxiliary variables. */
    char    *var;
    int      i, k;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */

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
    printf("1n====%d \n", n);

    b=(double*)malloc(n*sizeof(double));
    x=(double*)malloc(n*sizeof(double));
    ia=(int*)malloc((n+1)*sizeof(int));
    printf("2n====%d \n", n);
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%d", &row_index);
      ia[j] = row_index; 
    }
    fclose( f5);
    nnz = ia[n];
    printf("3nnz====%d \n", nnz);
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
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

    /* Set right hand side to i. */
    for (i = 0; i < n; i++) {
        b[i] = i;
    }

/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
    
    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
    if (error != 0) {
        printf("\nERROR in consistency of matrix: %d", error);
        exit(1);
    }

/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

    pardiso_chkvec (&n, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR  in right hand side: %d", error);
        exit(1);
    }

/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */

    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR right hand side: %d", error);
        exit(1);
    }
 
/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
    phase = 11; 

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
   
/* -------------------------------------------------------------------- */
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */    
    phase = 22;
    iparm[32] = 1; /* compute determinant */
    start = clock();
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("!!!!!!!!!!!!Numerical factorization took %f seconds to execute \n", cpu_time_used ); 
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */    
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
    phase = 33;

    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
   
    start = clock();
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("!!!!!!!!!!!!1Back substitution and iterative refinement took %f seconds to execute \n", cpu_time_used );
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }

    printf("\nSolve completed ... ");
    printf("\nThe solution of the system is: ");
    for (i = 0; i < n; i++) {
//        printf("\n x [%d] = % f", i, x[i] );
    }
    printf ("\n\n");


/* -------------------------------------------------------------------- */    
/* ... Inverse factorization.                                           */                                       
/* -------------------------------------------------------------------- */  

   if (solver == 0)
    {
    	printf("\nCompute Diagonal Elements of the inverse of A ... \n");
	phase = -22;
        iparm[35]  = 1; /*  no not overwrite internal factor L */ 

        start = clock();

        pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);

        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("!!!!!!!!!!!!1Inverse factorization took %f seconds to execute \n", cpu_time_used );
 /* print diagonal elements */
       for (k = 0; k < n; k++)
       {
            int j = ia[k]-1;      
  //          printf ("Diagonal element of A^{-1} = %d %d %32.24e\n", k, ja[j]-1, a[j]);
       }

    }   


/* -------------------------------------------------------------------- */    
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */ 
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */    
    phase = -1;                 /* Release internal memory. */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    if(ia!=NULL) free(ia);
    if(ja!=NULL) free(ja);
    if(a!=NULL) free(a);

    if(b!=NULL) free(b);
    if(x!=NULL) free(x);

    return 0;
} 
