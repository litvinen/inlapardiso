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
    fclose(fp);
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



/*convert graph representatio into CRS, 
 Note: graph model doesnt containdiagonal!*/
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


void convert2CSR_withdiag(psspmatrix S, graph_t * g, Qfunc_t Q, void *arg)
{
  int n = 0, nnz=0;
  int i, jj, j;
  int k = 0;
  
  printf("rows = %d, cols = %d\n", g->n, g->n);
  for (i = 0; i < g->n; i++) 
     nnz=nnz + g->nnbs[i];
  nnz = nnz + g->n;  // added number of diagonal elements, which are not in the graph
  S->nnz = nnz;
  S->n = g->n; 
  S->a = (double*)malloc(nnz * sizeof(double));
  S->ja =(int*)malloc(nnz*sizeof(int));
  S->ia=(int*)malloc(g->n*sizeof(int));
  S->ia[0]=0;
  
  for( k = 1; k <= g->n; k++)
  {
      S->ia[k] = S->ia[k-1] + g->nnbs[k-1]+1; 
  }
  for (i = 0; i < g->n; i++) {
    for (jj = 0; jj <= g->nnbs[i]; jj++)
    {    
        j = g->nbs[i][jj]-1;
        S->a[k] = Q(i, j, arg);
        k++;
    }
  }

  k=0;
  for (i = 0; i < g->n; i++) {
    for (j = 0; j <= g->nnbs[i]; j++)
    {    
        S->ja[k] = g->nbs[i][j]-1;
        k++;
    }
  }
  
}



void convert2CSR_alex(psspmatrix S, graph_t * g)
{
  int n = 0, nnz=0;
  int i, jj, j;
  int k = 0;
  double* values;
  
  printf("rows = %d, cols = %d\n", g->n, g->n);
  for (i = 0; i < g->n; i++) 
     nnz=nnz + g->nnbs[i];
  
  
  values=(double* )malloc(nnz*sizeof(double));
  for (i = 0; i < nnz; i++) 
    values[i] = 0.1*i;
    
  S->nnz = nnz;
  S->n = g->n; 
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
        S->ja[k] = g->nbs[i][j];
        k++;
    }
  }
  
}

/*convert graph representation from INLA into CRS and add a diagonal (contains in *diag), 
 Note: diagonal is added!*/
void convert2CSR_alex_diag(psspmatrix S, graph_t * g, double* diag)
{
  int n = 0, nnz=0;
  int i, jj, j;
  int k = 0;
  
  printf("rows = %d, cols = %d\n", g->n, g->n);
  for (i = 0; i < g->n; i++) 
     nnz=nnz + g->nnbs[i];
  
  nnz=nnz+g->n;
  
  S->nnz = nnz;
  S->n = g->n; 
  S->a = (double*)malloc(nnz * sizeof(double));
  S->ja =(int*)malloc(nnz*sizeof(int));
  S->ia=(int*)malloc(g->n*sizeof(int));
  S->ia[0]=0;
  
  for( k = 1; k <= g->n; k++)
     S->ia[k] = S->ia[k-1] + g->nnbs[k-1]+1;  //added 1 becouse of diagonal 
  for (i = 0; i < S->nnz; i++) 
     S->a[i] = diag[i];
  
  k=0;
  for (i = 0; i < g->n; i++) {
    for (j = 0; j <= g->nnbs[i]; j++)
    {    
        S->ja[k] = g->nbs[i][j];
        k++;
    }
  }
  
}


void test_alex_add_diag()
{
    graph_t *g;
    graph_t *newg;

    newg = (graph_t *) calloc(1, sizeof(graph_t));

    g = read_graph("minitest7.graph.txt");//working
    printf("Graph is readed succesfully! \n");
   // g = read_graph("minitest6.graph.txt"); //not working
    print_graph(g);
    /* Original matrix is
     0 1 2 0
     1 2 0 2
     2 3 0 1 3
     3 2 0 1
     */
    /* After adding main diagonal
     0 2 0 2 
     1 3 0 1 2
     2 4 0 1 2 3
     3 3 0 1 3
     */
    alex_add_diag(g, newg);
    print_graph(newg);

    
}


/*add diagonal to the INLA graph*/
void alex_add_diag(graph_t* g, graph_t* newg)
{
  /*Need to modify graph, so that is contains diagonal*/  
    int i, k, kk;
 //   graph_t *newg = (graph_t *) calloc(1, sizeof(graph_t));
 
    newg->n = g->n;
    
    newg->nnbs = (int *) calloc(newg->n, sizeof(int));
    newg->nbs = (int **) calloc(newg->n, sizeof(int *));

    for (k = 0; k < g->n; k++) 
       newg->nnbs[k] = g->nnbs[k] + 1;
    
    for (i = 0; i < g->n; i++) {
        newg->nbs[i] = (int *) calloc(newg->nnbs[i], sizeof(int));
        for (k = 0; k < g->nnbs[i]; k++) {
            if(g->nbs[i][k] < i)
              newg->nbs[i][k] = g->nbs[i][k];
            if(g->nbs[i][k] == i)
            {    
               printf("Attention, the innitial INLA graph cannot contain diagonal elements! \n");
            }   
            if(g->nbs[i][k] > i)
            {
              newg->nbs[i][k] = i;
              newg->nbs[i][k+1] = g->nbs[i][k];
            }   
        }
        if(g->nbs[i][g->nnbs[i]-1] < i)
           newg->nbs[i][newg->nnbs[i]-1] = i;
            
    }
    print_graph(newg);
}



void alex_add_diag_corrected(graph_t* g, graph_t* newg)
{
  /*Need to modify graph, so that is contains diagonal*/  
    int i, k, kk;
 //   graph_t *newg = (graph_t *) calloc(1, sizeof(graph_t));
 
    newg->n = g->n;
    
    newg->nnbs = (int *) calloc(newg->n, sizeof(int));
    newg->nbs = (int **) calloc(newg->n, sizeof(int *));

    for (k = 0; k < g->n; k++) 
       if(g->nnbs[k]==k) 
         newg->nnbs[k] = g->nnbs[k] + 1;
       else
         if(k<g->nnbs[k]) 
           newg->nnbs[k] = 1;
    
    for (i = 0; i < g->n; i++) {
        newg->nbs[i] = (int *) calloc(newg->nnbs[i], sizeof(int));
        for (k = 0; k < g->nnbs[i]; k++) {
            if(g->nbs[i][k] < i)
              newg->nbs[i][k] = g->nbs[i][k];
            if(g->nbs[i][k] == i)
            {    
               printf("Attention, the innitial INLA graph cannot contain diagonal elements! \n");
            }   
            if(g->nbs[i][k] > i)
            {
              newg->nbs[i][k] = i;
              newg->nbs[i][k+1] = g->nbs[i][k];
            }   
        }
        if(g->nbs[i][g->nnbs[i]-1] < i)
           newg->nbs[i][newg->nnbs[i]-1] = i;
            
    }
    print_graph(newg);
}


/*Print out arrays: a, ia, ja of a CSR matrix*/
void print_CSR(psspmatrix S)
{
    int k=0;
    printf("Array a[..]: \n"); 
    for( k = 0; k < S->nnz; k++)
       printf("%3.3g, ", S->a[k]); 
    printf(" \n"); 
    printf("Array ia[..]: \n"); 
    for( k = 0; k <= S->n; k++)
       printf("%d, ", S->ia[k]); 
    printf(" \n"); 
    printf("Array ja[..]: \n"); 
    for( k = 0; k < S->n; k++)
       printf("%d, ", S->ja[k]); 
    printf(" \n"); 
  
}


/* Allocation and initialization of PARDISO variables*/
int alex_mydata_initialization(pdata_storage mydata)
{
   int i=0;
   mydata->mtype=0;
   mydata->mnum=0;
   mydata->phase=0;
   mydata->idum=0;
   mydata->nrhs=0;
   mydata->msglvl=0;
   mydata->ddum=0.0;
   mydata->error=0;
   mydata->Q = (psspmatrix)malloc(sizeof(sspmatrix));
   mydata->Qinv = NULL;
   mydata->L = NULL;
//   mydata->Qinv = (psspmatrix)malloc(sizeof(sspmatrix));
//   mydata->L = (psspmatrix)malloc(sizeof(sspmatrix));
   mydata->mypermutation = NULL;
   
   for( i = 0; i < 64; i++)
   {
      mydata->iparm[i] = 0;
   }
   for( i = 0; i < 64; i++)
   {
      mydata->dparm[i] = 0.0; 
   }
}


/*Initialize Pardiso variables*/
int alex_initialization(int flag, pdata_storage mydata)
{
  int num_procs=0;
  char    *var;
  int i=0;
  
//flag==-1=="initial/common initialization of PARDISO"
//flag==0=="reordering"
//flag==1="symbolic factorization"
//flag==2="numerical factorization"
//flag==3="Cholesky factorization"
//flag==4="Partial inversion"
   printf("flag= %d \n", flag);

   

   
   if(flag==-1) //general initialization
    { 
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
    }
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


     pardiso_chkmatrix  (&(mydata->mtype), &(mydata->Q->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->error));
     if (mydata->error != 0) {
        printf("\nERROR in consistency of matrix: %d", mydata->error);
        exit(1);
     }
   

   if(flag==0) //0="reordering"
   {
     /* -------------------------------------------------------------------- */
     /* ..  Setup Pardiso control parameters.                                */
     /* -------------------------------------------------------------------- */
      mydata->phase = 11; // analysis
      mydata->error = 0;
      mydata->iparm[4] = 0; //set to 1 if you want to provide your own permutation array
    

   }
   if(flag==1) //1="symbolic factorization"
   {
      mydata->phase = 11; // analysis
   }
   if(flag==2) //2="numerical factorization"
   {
      mydata->phase = 22; //Numerical factorization
      mydata->iparm[32] = 1; /* compute determinant */
   }
   if(flag==3) //3="Cholesky factorization"
   {
      mydata->L = (psspmatrix)malloc(sizeof(sspmatrix));

      mydata->phase = 12; // Analysis, numerical factorization
      //  IPARM (26) — Splitting of Forward/Backward Solve.
      
   }
   if(flag==4) //4="Selected inversion"
   {
      if(mydata->Qinv!=NULL)
          printf("Possible error! mydata->Qinv!=NULL. Memory will be allocated the second time!");
      mydata->Qinv = (psspmatrix)malloc(sizeof(sspmatrix));
      mydata->phase = -22; /*  do not overwrite internal factor L with selected inversion*/ 
      mydata->iparm[35]  = 1; /*  do not overwrite internal factor L with selected inversion*/ 
/*It will instead allocate additional memory of size of the numbers of
elements in L and U to store these inverse elements. 
*/
//  VERY DANGEROUS OPTION    mydata->iparm[36]  = 1; /* PARDISO will return the selected inverse elements
//A^−1_ij in full symmetric triangular CSR format.*/ 
   }
   if(flag==5) //4="Solve for x"
   {
      mydata->iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
      mydata->phase = 33;  //solve
   
   }
   if(flag==6) //6="Cholesky factorization"
   {
      mydata->phase = 12; // Analysis, numerical factorization
      //  IPARM (26) — Splitting of Forward/Backward Solve.
      
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

/* Read sparse matrix in CSR format */
void read_CSR_matrix( psspmatrix Q, char *filename)
{
    int    n = 0;
    int    nnz = 0;
    FILE *f5;
    int anzv=0, j=0, row_index=0, col_index=0;
    double elem=0.0;
    double* a=NULL;
    int* ia=NULL;
    int* ja=NULL;
    
    f5 = fopen( filename, "r");
    if( f5 == NULL)
    {
      printf( "Can't open the file\n");
      exit (1);
    }
    fscanf( f5, "%d\n", &anzv);
    n = anzv-1;
    //printf("n=%d \n", anzv);
    Q->n=n;
    Q->ia=(int*)malloc((n+1)*sizeof(int));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%d", &row_index);
      //printf("row_index=%d \n", row_index);
      Q->ia[j] = row_index; 
    }
    
    fscanf( f5, "%d\n", &anzv); 
    Q->nnz=anzv;
    Q->ja=(int*)malloc((anzv+1)*sizeof(int));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%d", &col_index);
      Q->ja[j] = col_index; 
    }
    fscanf( f5, "%d\n", &anzv); 
    //printf("anzv=%d \n", anzv);
    Q->a=(double*)malloc((anzv+1)*sizeof(double));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%lf", &elem);
      Q->a[j] = elem; 
    }
    //Q->ia = ia;
    //Q->ja = ja;
    //Q->a  = a;
    
    fclose( f5);
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
    printf("n=%d \n", anzv);
    Q->n=n;
    Q->ia=(int*)malloc((n+1)*sizeof(int));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%d", &row_index);
      //printf("row_index=%d \n", row_index);
      Q->ia[j] = row_index; 
    }
    fclose( f5);
    
    f5 = fopen( "ja.txt", "r");
    if( f5 == NULL)
    {
      printf( "Can't open the file\n");
      exit (1);
    }
    fscanf( f5, "%d\n", &anzv); 
    Q->nnz=anzv;
    Q->ja=(int*)malloc((anzv+1)*sizeof(int));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%d", &col_index);
      Q->ja[j] = col_index; 
    }
    fclose( f5);
    f5 = fopen( "a.txt", "r");
    if( f5 == NULL)
    {
      printf( "Can't open the file\n");
      exit (1);
    }
    fscanf( f5, "%d\n", &anzv); 
    printf("anzv=%d \n", anzv);
    Q->a=(double*)malloc((anzv+1)*sizeof(double));
    for( j = 0; j < anzv; j++)
    {
      fscanf( f5, "%lf", &elem);
      Q->a[j] = elem; 
    }
    //Q->ia = ia;
    //Q->ja = ja;
    //Q->a  = a;
    
    fclose( f5);
}   





// reordering and symbolic fact., computed just once
int alex_reordering(pdata_storage mydata, psgraph_t mgraph)
{
    clock_t start, end;
    double cpu_time_used;
//IPARM (5) — User permutation.
    int k=0;
    
    for(k=0; k < mydata->Q->n; k++)
        printf("Before mypermutation[%d]=%d\n", k, mydata->mypermutation[k]);
    
    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
	     &(mydata->Q->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, mydata->mypermutation, &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error), mydata->dparm);

    for(k=0; k < mydata->Q->n; k++)
        printf("After mypermutation[%d]=%d\n", k, mydata->mypermutation[k]);
    
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
   
    /*pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
	     &(mydata->Q->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error), mydata->dparm);
    
    if (mydata->error != 0) {
        printf("\nERROR during symbolic factorization: %d", mydata->error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", mydata->iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", mydata->iparm[18]);
*/
  return 0;
}


/*Create a full copy of a CSR matrix*/
void alex_CSRmatrix_copy(psspmatrix  source, psspmatrix  destin, pdata_storage mydata)
{
    int i=0;
    pardiso_chkmatrix(&(mydata->mtype), &(source->n), source->a,  source->ia,  source->ja, &( mydata->error));
    if ( mydata->error!= 0) {
        printf("\n source  ERROR in consistency of matrix: %d",  mydata->error);
        exit(1);
    }
    else
      printf("\n source is fine! \n");        

    
   destin->n = source->n ;
   destin->nia = source->nia ;
   destin->nnz = source->nnz ;
   destin->logdet = source->logdet ;
   destin->ia=(int*)malloc((1+source->n+1) * sizeof(int)); 
   destin->ja=(int*)malloc(source->nnz * sizeof(int)); 
   destin->a=(double*)malloc(source->nnz * sizeof(double)); 
   //memcpy(destin->ia, source->ia, (source->n+1)*sizeof(int) );
   //memcpy(destin->ja, source->ja, source->nnz*sizeof(int) );
   //memcpy(destin->a, source->a, source->nnz*sizeof(double) );
   for (i = 0; i <= source->n; i++) 
       destin->ia[i] = source->ia[i];
   for (i = 0; i < source->nnz; i++) {
       destin->ja[i] = source->ja[i];
       destin->a[i] = source->a[i];
   }    
   
   
    pardiso_chkmatrix(&(mydata->mtype), &(destin->n), destin->a,  destin->ia,  destin->ja, &( mydata->error));
    if ( mydata->error!= 0) {
        printf("\n !!!1111ERROR in consistency of matrix: %d",  mydata->error);
        exit(1);
    }
    else
      printf("\n destin is fine! \n");        
        

}

//Numerical Cholesky factorization
//int alex_chol(pdata_storage mydata, pINLA_mtx Q)
void alex_chol(pdata_storage mydata)
{
   printf("Copy matrix Q into matrix L. MAY BE PARDISO DO IT AUTOMATICALLY !?!?!? \n");
   alex_CSRmatrix_copy(mydata->Q, mydata->L,  mydata);
   pardiso_chkmatrix  (&(mydata->mtype), &(mydata->L->n), mydata->L->a, mydata->L->ia, mydata->L->ja, &(mydata->error));
    if (mydata->error != 0) {
        printf("\n  !!!ERROR in consistency of matrix: %d", mydata->error);
        exit(1);
    }
    else
      printf(" L is fine\n");
    
    
    
    printf("Copy is done. \n");
    mydata->iparm[32] = 1; /* compute determinant */
    //    printf("Matrix L is not empty !!! \n");
    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->L->n), mydata->L->a, mydata->L->ia, mydata->L->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error),  mydata->dparm);
  
    if (mydata->error != 0) {
        printf("\nERROR during numerical factorization: %d", mydata->error);
        exit(2);
    }
    printf("\n Cholesky factorization completed ...\n ");
    printf("\n !!! Determinant = %3.3g ...\n ", mydata->dparm[32]);
    

}




//Solve linear system Lx=b
void alex_forward_subst(pdata_storage mydata,  double* rhs,  double* y)
{
/* ..  Back substitution and iterative refinement.                      */
   pardiso ( mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->L->n), mydata->L->a, mydata->L->ia, mydata->L->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), rhs, y, &(mydata->error),  mydata->dparm);
   if (mydata->error != 0) {
        printf("\nERROR during solution: %d", mydata->error);
        exit(3);
    }
    printf("\nSolve Ly=rhs completed ... ");
    printf ("\n\n");

}


/* -------------------------------------------------------------------- */    
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
//Solve linear system L^Tx=b
void alex_back_subst(pdata_storage mydata,  double* y,  double* x)
{
    clock_t start, end;
    double cpu_time_used;
    int i=0;
//On entry: Solve a system A^T x = b by using the factorization of A.
//The default value of IPARM(12) is 0. IPARM(12)=1 indicates that you would like to solve
//transposed system A^T x = b.
   mydata->iparm[11]=1;

   /* start = clock();*/
    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->L->n), mydata->L->a, mydata->L->ia, mydata->L->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), y, x, &(mydata->error),  mydata->dparm);
    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("!!!!!!!!!!!!1Back substitution and iterative refinement took %f seconds to execute \n", cpu_time_used );
    if (mydata->error != 0) {
        printf("\nERROR during solution: %d", mydata->error);
        exit(3);
    }

    printf("\nSolve L^Tx=y is also completed ... ");
    printf("\nThe solution of the system is: ");
    for (i = 0; i < mydata->Q->n; i++) {
        printf("\n x [%d] = % f", i, x[i] );
    }
    printf ("\n\n");
}


//Solve linear system LL^Tx=b, return x
void alex_solve_LLTx(pdata_storage mydata,  double* rhs,  double* x)
{
  double* y=NULL;
  int j=0, k=0;
  y = (double*)malloc( mydata->Q->n*sizeof(double)) ;
  //  IPARM (26) — Splitting of Forward/Backward Solve.
  mydata->iparm[25] = 1; // IPARM(25) = 1 indicates that a forward solve step is performed with the factors L
  alex_forward_subst(mydata,  rhs, y); // L y = rhs
  mydata->iparm[25] = 2; // IPARM(25) = 2 indicates that a backward solve step is performed with the factors U or L T
  alex_back_subst(mydata, y, x); // L^T x = y
  
  free(y);

}

//Computing partial inverse, means solving linear system foer specific RHS
void alex_inv(pdata_storage mydata)
{
  psspmatrix Q = mydata->Q;
  
  /*On entry: IPARM(36) will control the selected inversion process based on the internal L and
  U factors. If IPARM(36) = 0 and PHASE = -22, PARDISO will overwrite these factors with
  the selected inverse elements of A −1 . If IPARM(36) = 1 and PHASE = -22, PARDISO will
  not overwrite these factors. It will instead allocate additional memory of size of the numbers of
  elements in L and U to store these inverse elements.*/
  int k=0, i=0, j=0;
  double* x;
  double* b;
  x = (double*)malloc( mydata->Q->n*sizeof(double)) ;
  b = (double*)malloc( mydata->Q->n*sizeof(double)) ;
  alex_CSRmatrix_copy(mydata->Q, mydata->Qinv,  mydata);
  
   for (i = 0; i < mydata->Q->n; i++) {
      b[i] = i;
      x[i] = 0.0;
   }
   if (mydata->solver == 0)
   {
        mydata->phase = -22;
        mydata->iparm[35]  = 1; /*  do no not overwrite internal factor L */ 
        pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase), &(mydata->Qinv->n), mydata->Qinv->a, mydata->Qinv->ia, mydata->Qinv->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), b, x, &(mydata->error),  mydata->dparm);

       //        printf("!!!!!!!!!!!!1Inverse factorization took %f seconds to execute \n", cpu_time_used );
       printf ("Diagonal element of A^{-1}\n");
       for (k = 0; k <7; k++)
       {
            j = mydata->Q->ia[k]-1;
            printf ("A^{-1}(%d, %d)= %32.24e\n", k, mydata->Qinv->ja[j]-1, mydata->Qinv->a[j]);
       }
   } 
  
}

//Compute log-det
double alex_log_det(pdata_storage mydata)
{
   double logdet=0.0;

   pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->Q->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error),  mydata->dparm);
   printf("\n !!!NEW!!! Determinant = %3.3g ...\n ", mydata->dparm[32]);

/*
  IPARM (33) — Determinant of a matrix.
  Input
  On entry: IPARM(33)=1 will compute the determinant of matrices and will return in DPARM(33)
  the real part of the determinant and, if necessary, in DPARM(32) the complex part of the deter-
  minant.
  On output: The parameter returns the natural logarithm of the determinant of a sparse matrix
  A. The default value of IPARM(33) is 0.
*/

  mydata->Q->logdet=   mydata->dparm[32];
  return mydata->dparm[32];
}


//Finalize and deallocate memory

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */    
void alex_finalize(pdata_storage mydata)
{
    int phase = -1;                 /* Release internal memory. */
    
    if(mydata!=NULL)
    {
  /*      if(mydata->Q!=NULL)
            alex_clean_CRS_matrix(mydata->Q);*/
        if(mydata->Qinv!=NULL)
            alex_clean_CRS_matrix(mydata->Qinv);
        if(mydata->L!=NULL)
           alex_clean_CRS_matrix(mydata->L);
        if(mydata->mypermutation!=NULL)
        {    
          free(mydata->mypermutation);
          mydata->mypermutation=NULL;
        }
    }
    pardiso (mydata->pt, &(mydata->maxfct), &(mydata->mnum), &(mydata->mtype), &(mydata->phase),
             &(mydata->Q->n), &(mydata->ddum), mydata->Q->ia, mydata->Q->ja, &(mydata->idum), &(mydata->nrhs),
             mydata->iparm, &(mydata->msglvl), &(mydata->ddum), &(mydata->ddum), &(mydata->error),  mydata->dparm);


}

void  alex_clean_CRS_matrix(psspmatrix A)
{
    if(A!=NULL)
    {
        if(A->a!=NULL)
            free(A->a);
        if(A->ia!=NULL)
            free(A->ia);
        if(A->ja!=NULL)
            free(A->ja);
        free(A);
        A=NULL;

    }
}

int alex_clean_mydata(int flag, pdata_storage  mydata)
{
   /* To clean PARDISO flags after each step... */
//flag==-1=="initial/common initialization of PARDISO"
//flag==0=="reordering"
//flag==1="symbolic factorization"
//flag==2="numerical factorization"
//flag==3="Cholesky factorization"
//flag==4="Partial inversion"
//4="Solve for x"
//6="Cholesky factorization"
   if(flag==0) //0="reordering"
   {
     /* -------------------------------------------------------------------- */
     /* ..  Setup Pardiso control parameters.                                */
     /* -------------------------------------------------------------------- */
      mydata->phase = 0; // analysis
      mydata->error = 0;
      mydata->iparm[4] = 0; //set to 1 if you want to provide your own permutation array
   }
   if(flag==1) //1="symbolic factorization"
   {
      mydata->phase = 0; // analysis
   }
   if(flag==2) //2="numerical factorization"
   {
      mydata->phase = 0; //Numerical factorization
      mydata->iparm[32] = 0; /* compute determinant */
   }
   if(flag==3) //3="Cholesky factorization"
   {
      mydata->phase = 0; // Analysis, numerical factorization
      //  IPARM (26) — Splitting of Forward/Backward Solve.
      
   }
   if(flag==4) //4="Selected inversion"
   {
      mydata->phase = 0; /*  do not overwrite internal factor L with selected inversion*/ 
      mydata->iparm[35]  = 0; /*  do not overwrite internal factor L with selected inversion*/ 
/*It will instead allocate additional memory of size of the numbers of
elements in L and U to store these inverse elements. 
*/
//  VERY DANGEROUS OPTION    mydata->iparm[36]  = 1; /* PARDISO will return the selected inverse elements
//A^−1_ij in full symmetric triangular CSR format.*/ 
   }
   if(flag==5) //4="Solve for x"
   {
      mydata->iparm[7] = 0;       /* Max numbers of iterative refinement steps. */
      mydata->phase = 0;  //solve
   
   }
   if(flag==6) //6="Cholesky factorization"
   {
      mydata->phase = 0; // Analysis, numerical factorization
      //  IPARM (26) — Splitting of Forward/Backward Solve.
      
   }
   
}





int test_conversion()
{
    graph_t *g;
    psspmatrix S;
  
    S=(psspmatrix )malloc(sizeof(sspmatrix ));
    //g = read_graph("germany.graph.txt");
    //g = read_graph("minitest4.graph.txt");
    //g = read_graph("minitest2.graph.txt");
    g = read_graph("minitest3.graph.txt");
    /* CRS matrix is
     1 0 0 0
     5 8 0 0
     0 0 3 0
     0 6 0 1
     */
    print_graph(g);
   // print_Q(g, Qdemo, (void *) g);
    //convert2CSR(S, g, Qdemo, (void *) g);

    convert2CSR_alex(S, g);
    print_CSR(S);
    /*
     should obtain
     A= [5 8 3 6]
     IA=[0 0 2 3 4]
     JA=[0 1 2 1]
     */
    
    return (0);
}


void store_CSR_matrix(psspmatrix S, char *filename)
{
    int k=0;
    FILE *fp = NULL;
    
    fp = fopen(filename, "w");
    
    fprintf(fp, "%d\n", S->nnz);
    fprintf(fp, "%d\n", S->n);
    
    for( k = 0; k < S->nnz; k++)
       fprintf(fp, "%3.3g\n ", S->a[k]); 
    fprintf(fp, " \n"); 
    //fprintf(fp, "Array ia[..]: \n"); 
    for( k = 0; k <= S->n; k++)
       fprintf(fp, "%d\n", S->ia[k]); 
    fprintf(fp, " \n"); 
    //fprintf(fp, "Array ja[..]: \n"); 
    for( k = 0; k < S->nnz; k++)
       fprintf(fp, "%d\n", S->ja[k]); 
    fclose(fp);
  
}


/*Restore CSR matrix from file, assuming that memory for a, ia, ja is NOT allocated*/
void restore_CSR_matrix(psspmatrix S, char *filename)
{
    int i=0, ip=0;
    double dp=0.0;
    FILE *fp = NULL;
    
    fp = fopen(filename, "r");
    
    fscanf( fp, "%d\n", &(S->nnz));
    fscanf( fp, "%d\n", &(S->n));
    
    S->ia = (int*)malloc((S->n+1)*sizeof(int));
    S->ja = (int*)malloc(S->nnz*sizeof(int));
    S->a  = (double*)malloc(S->nnz*sizeof(double));

    for( i = 0; i < S->nnz; i++)
    {
      fscanf( fp, "%lf\n ", &dp);
      S->a[i] = dp;
    }
    fscanf( fp, "\n");
    for( i = 0; i <= S->n; i++)
    {
      fscanf( fp, "%d\n ", &ip);
      S->ia[i] = ip;
    }
    fscanf( fp, "\n");
    for( i = 0; i < S->nnz; i++)
    {
      fscanf( fp, "%d\n ", &ip);
      S->ja[i] = ip;
    }

    fclose(fp);
  
}


/*Restore pardiso internal structures from a file*/
void alex_pardiso_restore(pdata_storage mydata, char *filename)
{
    int i=0, ip=0;
    double dp=0.0;
    FILE *fp = NULL;
    
    fp = fopen(filename, "r");
    fscanf( fp, "%d\n", &mydata->mtype);
    fscanf( fp, "%d\n", &mydata->mnum);
    fscanf( fp, "%d\n", &mydata->phase);
    fscanf( fp, "%d\n", &mydata->idum);
    fscanf( fp, "%d\n", &mydata->nrhs);
    fscanf( fp, "%d\n", &mydata->msglvl);
    fscanf( fp, "%lf\n", &mydata->ddum);
    fscanf( fp, "%d\n", &mydata->error);
    for( i = 0; i < 64; i++)
    {
      fscanf( fp, "%d\n ", &ip);
      mydata->iparm[i] = ip;
    }
    for( i = 0; i < 64; i++)
    {
      fscanf( fp, "%lf\n ", &dp);
      mydata->dparm[i] = dp; 
    }
    for( i = 0; i < mydata->Q->n; i++)
    {
      fscanf( fp, "%d\n ", &ip);
      mydata->mypermutation[i] = ip;
    }
    fclose(fp);
    //filename = "matrixCSR.txt";
    
    /* Read sparse matrix in matlab format */
    if(mydata->Q!= NULL)
       mydata->Q=(psspmatrix )malloc(sizeof(sspmatrix ));
    restore_CSR_matrix( mydata->Q,  "alex_Q.txt");
    store_CSR_matrix(mydata->Q, "alex_Q2.txt");
    
    //if(mydata->Qinv!= NULL)
    //   mydata->Qinv=(psspmatrix )malloc(sizeof(sspmatrix ));
    //restore_CSR_matrix( mydata->Qinv,  "alex_Qinv.txt");
    //store_CSR_matrix(mydata->Qinv, "alex_Qinv.txt");
    //if(mydata->L!= NULL)
    //   mydata->L=(psspmatrix )malloc(sizeof(sspmatrix ));
    //restore_CSR_matrix( mydata->L,  "alex_L.txt");
    
    
    
}

/*Store internal structures from pardiso to a file*/
void alex_pardiso_store(pdata_storage mydata, char *filename)
{
    FILE *fp = NULL;
    int i=0;
    fp = fopen(filename, "w");
    
    /*
    mydata->mtype=1;
    mydata->mnum=2;
    mydata->phase=3;
    mydata->idum=4;
    mydata->nrhs=5;
    mydata->msglvl=6;
    mydata->ddum=7.7;
    mydata->error=8;
    */
    
    fprintf(fp, "%d\n",  mydata->mtype);
    fprintf(fp, "%d\n",  mydata->mnum);
    fprintf(fp, "%d\n",  mydata->phase);
    fprintf(fp, "%d\n",  mydata->idum);
    fprintf(fp, "%d\n",  mydata->nrhs);
    fprintf(fp, "%d\n",  mydata->msglvl);
    fprintf(fp, "%12.12g\n",  mydata->ddum);
    fprintf(fp, "%d\n",  mydata->error);
    for(i=0; i<64; i++)
       fprintf(fp, "%d\n", mydata->iparm[i]);
    for(i=0; i<64; i++)
       fprintf(fp, "%12.12g\n", mydata->dparm[i]);
    for(i=0; i<mydata->Q->n; i++)
       fprintf(fp, "%d\n", mydata->mypermutation[i]);
    fclose(fp);
    
        
        
    //filename = "matrixCSR.txt";
    if(mydata->Q!= NULL)
        store_CSR_matrix(mydata->Q, "alex_Q.txt");
    if(mydata->Qinv!= NULL)
        store_CSR_matrix(mydata->Qinv, "alex_Qinv.txt");
    if(mydata->L!= NULL)
        store_CSR_matrix(mydata->L, "alex_L.txt");
    
}

void alex_test_compare_with_pardiso(pdata_storage mydata)
{
       /* Matrix data. */
    int    n = 8;
    int* ia=NULL;
    int* ja=NULL;
    double* a=NULL;
    int* mypermutation;
    
    mypermutation = (int*)malloc(9*sizeof(int));
    ia = (int*)malloc(9*sizeof(int));
    ja = (int*)malloc(18*sizeof(int));
    a  = (double*)malloc(18*sizeof(double));
    
    ia[ 0] = 0;
    ia[ 1] = 4;
    ia[ 2] = 7;
    ia[ 3] = 9;
    ia[ 4] = 11;
    ia[ 5] = 14;
    ia[ 6] = 16;
    ia[ 7] = 17;
    ia[ 8] = 18;

//    ia[ 9] = { 0, 4, 7, 9, 11, 14, 16, 17, 18 };
    ja[0] =  0;
    ja[1] =  2;
    ja[2] =  5;
    ja[3] =  6;
    ja[4] =  1;
    ja[5] =  2;
    ja[6] =  4;
    ja[7] =  2;
    ja[8] =  7;
    ja[9] =  3;
    ja[10] =  6;
    ja[11] =  4;
    ja[12] =  5;
    ja[13] =  6;
    ja[14] =  5;
    ja[15] =  7;
    ja[16] =  6;
    ja[17] =  7;
    
    a[0] = 7.0;
    a[1] = 1.0;
    a[2] = 2.0;
    a[3] = 7.0;
    a[4] = -4.0;
    a[5] = 8.0;
    a[6] = 2.0;
    a[7] = 1.0;
    a[8] = 5.0;
    a[9] = 7.0;
    a[10] = 9.0;
    a[11] = 5.0;
    a[12] = 1.0;
    a[13] = 5.0;
    a[14] = 0.0;
    a[15] = 5.0;
    a[16] = 11.0;
    a[17] = 5.0;

    int     nnz = ia[n];
    int     mtype = -2;        /* Real symmetric matrix */
//    Q=(psspmatrix )malloc(sizeof(sspmatrix ));
    
    mydata->mypermutation=mypermutation;
    mydata->Q->nnz=nnz;
    mydata->Q->n=n;
    mydata->Q->ia=ia;
    mydata->Q->a=a;
    mydata->Q->ja=ja;
    mydata->mtype=mtype;
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
    void    *pt[64]; /* Internal solver memory pointer pt,                  */
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
    graph_t *g;
    
    
    printf("!!!TEST!!!!\n");
    test_alex_add_diag();
    printf("!!!TEST ENDED!!!!\n");
    
    
    printf("!!!Note that in C programming : iparm[3]=number_of_processors from manual is iparm[2] in C program!!!!\n");
    mydata= (pdata_storage)malloc(sizeof(sdata_storage ));
    alex_mydata_initialization(mydata);

   // test_conversion();

    solver=0;/* use sparse direct solver */
    mtype=-2; // real symmetric 
    nrhs=1;
  
    
    mydata->mtype = mtype;
    mydata->nrhs = nrhs;
    mydata->solver = solver;
    pardisoinit(mydata->pt,  &(mydata->mtype), &(mydata->solver), mydata->iparm, mydata->dparm,  &(mydata->error)); 
   // One of the options is just to read CRS matrix from a file(s) read_sparse_matrix(mydata->Q);
   //Another option is to take it from INLA
    
    
    //g = read_graph("germany.graph.txt");
    //g = read_graph("minitest3.graph.txt");
    //print_graph(g);
   // convert2CSR_alex(mydata->Q, g);
   //mydata->Q = (psspmatrix)malloc(sizeof(sspmatrix));
   //mydata->Qinv = (psspmatrix)malloc(sizeof(sspmatrix));
   //mydata->L = (psspmatrix)malloc(sizeof(sspmatrix));
   alex_test_compare_with_pardiso(mydata); //Just to test .

    //print_CSR(mydata->Q);

    // print_Q(g, Qdemo, (void *) g);
    //convert2CSR_alex_diag(mydata->Q, g, values);

   //alex_pardiso_store(mydata, "a.txt");
    
    
    x = (double*)malloc( mydata->Q->n*sizeof(double)) ;
    b = (double*)malloc( mydata->Q->n*sizeof(double)) ;
    for( j = 0; j < mydata->Q->n; j++)
    {
      mydata->mypermutation[j] = j+1;
      x[j]=0.0; 
      b[j]=j; 
    } 
       
    convert_C2F(mydata->Q);
    pardiso_chkmatrix  (&(mydata->mtype), &(mydata->Q->n), mydata->Q->a, mydata->Q->ia, mydata->Q->ja, &(mydata->error));
    if (mydata->error != 0) {
        printf("\nERROR in consistency of matrix: %d", mydata->error);
        exit(1);
    }
    
    alex_initialization(0, mydata); //flag==0=="reordering"
    alex_reordering(mydata, mgraph);
    alex_clean_mydata(0, mydata);

    alex_initialization(1, mydata); //flag==1="symbolic factorization"
    alex_initialization(2, mydata); //flag==2="numerical factorization"
    alex_symbolic_factorization(mydata);
    alex_clean_mydata(1, mydata);
    alex_clean_mydata(2, mydata);

    alex_initialization(3, mydata); //flag==3="Cholesky factorization"
    alex_chol(mydata);  //also computes log determinant
    alex_clean_mydata(3, mydata);

    alex_initialization(6, mydata); //flag==6="Compute only log-determinant (Cholesky is also computed)"
    alex_log_det(mydata);
    alex_clean_mydata(6, mydata);

    alex_initialization(5, mydata); //flag==5="Solution"
    alex_solve_LLTx(mydata, b, x);
    alex_clean_mydata(5, mydata);
    alex_pardiso_store(mydata, "alex_test_output.txt");
    alex_pardiso_restore(mydata, "alex_test_output.txt");

    alex_initialization(4, mydata); //flag==4="Partial inversion"
    alex_inv(mydata);
    convert_F2C(mydata->Q);
    alex_clean_mydata(4, mydata);
    alex_finalize(mydata);
    if(x!=NULL) free(x);
    if(b!=NULL) free(b);
    printf("alex_clean_mydata: Memory is deallocated.\n");
    
    
    return 0;
}


