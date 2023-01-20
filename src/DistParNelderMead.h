/*
 * DistParNelderMead.hpp
 *
 * Implements MPI based distributed memory parallel NelderMead simplex method.
 *
 *  Created on: May 10, 2011
 *      Author: kyleklein
 *  Modified on Jan 20, 2023 for R
 *      Author: chrischeungnf
 */

#ifndef NELDERMEAD_HPP_
#define NELDERMEAD_HPP_
#define SIMPLEX(i,j) simplex[((indices[(i)])*dimension) + (j) ]
#define ALPHA (1.0)
#define BETA (.5)
#define GAMMA (1.0)
#define TAU (.5)
#include <Rcpp.h>
#include <mpi.h>
#include <cstring> // memmove
#include <cmath> // fabs
#include <cfloat> // DBL_MAX
#include <algorithm> // min
#define DEBUG 0
#if DEBUG
  #include <fstream> // For logging
  #include <string> // to_string
#endif

class DistParNelderMead {
public:
	/**
	 * Given initial guess and a pointer to an objective function.
	 *
	 * simplex: Array of doubles size dimension*(dimension+1)
	 * dimensions: Dimension of each of the dimension+1 vectors.
	 * obj_function: Pointer to the objective function, takes as argument
	 * 	            a vector and an argument list, should return a double.
	 */
	DistParNelderMead(Rcpp::NumericVector& guess, Rcpp::Function* obj_function, Rcpp::List& obj_function_arg,
                   int size, int rank);
  
	/*
	 * Deletes user passed simplex as well as all allocated memory.
	 */
	~DistParNelderMead();
	/**
	 * Find the point which minimizes the objective function, and return
	 * an array of dimension doubles. User is responsible to free that memory.
	 *
	 * Will return answer if no improvement for 10 consecutive iterations, or if max_iterations > 0, then after
	 * max_iterations, whichever comes first.
	 */
	Rcpp::List solve(int max_iterations);
	// //Set alpha, otherwise assumed to be ALPHA
	// void setAlpha(double alpha);
	// //Set beta, otherwise assumed BETA
	// void setBeta(double beta);
	// //Set gamma, otherwise assumed GAMMA
	// void setGamma(double gamma);
	// //Set tau, otherwise assumed to be TAU
	// void setTau(double tau);
	// //Set minimimum improvement to do restart after some number of iterations
	// void setRestartCriterion(int iterations, double improvement);
	
private:
  Rcpp::NumericVector C_Array_to_R_Vec(double *vector, int dimension);
  double point_value(double* point);
  void create_simplex(double *guess);
	void centroid();
	void reflection();
	void expansion();
	void contraction();
	
	void global_best(MPI_Comm comm);
	void minimize();
	void daxpy(double *result, double scalar1, double *a, double scalar2,
			double *b, int length);
	void sort_simplex();
	void update(double *vector, int index);
	
	double *simplex, *M, *AR, *AE, *AC, *global_best_par;
	double *obj_function_results;
	double alpha, beta, gamma, tau, fAR, fAE, fAC, best, prev_best;
	int *indices;
	int dimension, points_on_proc;
	int rank, size_used, points_per_iter, current_point;
	int updated, convergence;
	MPI_Comm work_comm;
	Rcpp::Function* obj_function;
	Rcpp::List obj_function_arg;

  #if DEBUG
  	void print_simplex();
  	void print_simplex_with_value();
	  std::ofstream logfile;
  #endif
};

class IndexSorter {
public:
	IndexSorter(double *arg) :
			obj_function_results(arg) {};

	bool operator()(int i, int j) {
		return obj_function_results[i] < obj_function_results[j];
	}
private:
	double *obj_function_results;
};

#endif /* NELDERMEAD_HPP_ */
