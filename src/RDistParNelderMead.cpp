/*
 * RDistParNelderMead.cpp
 *
 * Call MPI based distributed memory parallel NelderMead simplex C++ method from R.
 *
 *  Created on: Jan 20, 2023 for R
 *      Author: chrischeungnf
 */
#include "DistParNelderMead.h"
using namespace Rcpp;

// [[Rcpp::export]]
List RDistParNelderMead(NumericVector par, Function fn, List fn_arg, int max_iterations=-1) {
  int mpi_inited;
  // Initialize the MPI environment
  MPI_Initialized(&mpi_inited);
  if(!mpi_inited) {
    MPI_Init(NULL, NULL);
  }
  
  int size, rank;
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  DistParNelderMead *solver = new DistParNelderMead(par, &fn, fn_arg, size, rank);
  
  //double t1, t2;
  //t1 = MPI_Wtime();
  List res = solver->solve(max_iterations);
  //t2 = MPI_Wtime();
  //if (rank == 0)
  //  Rcout << "Elapsed time during solve: " << t2 - t1 << "\n";
  
  delete solver;
  return res;
}

// [[Rcpp::export]]
void RMPI_Finalize() {
  MPI_Finalize();
}
