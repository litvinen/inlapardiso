// This file contains interface between sparse solver Pardiso and INLA library
// 14 January 2018
// KAUST, Bayesian Computational statistics and modeling, CEMSE, KAUST, 
// https://bayescomp.kaust.edu.sa/Pages/Home.aspx
// Dr. Alexander Litvinenko and Prof. Haavrd Rue


typedef struct data_storage sdata_storage;
typedef sdata_storage* pdata_storage;

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


/*Coincidence matrix (GRAPH)*/
typedef struct coin_matrix scoin_matrix;
typedef scoin_matrix* pcoin_matrix;

/*The structure which contains all required data*/
struct coin_matrix{
  int rows;
  int cols;
  int* ia;
  int* ja;
};



// reordering, computed just once
int alex_reordering(pdata_storage mydata_storage, pgraph* mgraph)
{
  return 0;
}

//compute the symbolic factorization, computed just once
int alex_symbolic_factorisation(data_storage* mydata_storage)
{
  return 0;
}

//Cholesky factorization
int alex_chol(data_storage* mydata_storage, pmatrix Q)
{
  return 0;
}

//Solve linear syatem Lx=b
int alex_solve_Lx(data_storage* mydata_storage,  double* b)
{
  return 0;
}


//Solve linear system L^Tx=b
int alex_solve_LTx(data_storage* mydata_storage,  double* b)
{
  return 0;
}


//Solve linear system LL^Tx=b
int alex_solve_LLTx(data_storage* mydata_storage,  double* b)
{
  return 0;
}

//Computing partial inverse, means solving linear system foer specific RHS
int alex_inv(data_storage* mydata_storage,  pmatrix Q)
{
  return 0;
}

//Compute log-det
int alex_log_det(data_storage* mydata_storage)
{
  return 0;
}





