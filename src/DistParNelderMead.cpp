/*
 * DistParNelderMead.cpp
 *
 * Implements MPI based distributed memory parallel NelderMead simplex method.
 *
 *  Created on: May 10, 2011
 *      Author: kyleklein
 *  Modified on Jan 20, 2023 for R
 *      Author: chrischeungnf
 */

#include "DistParNelderMead.h"

DistParNelderMead::DistParNelderMead(Rcpp::NumericVector& guess, Rcpp::Function* obj_function, Rcpp::List& obj_function_arg,
                                     int size, int rank) {
  dimension = guess.size();
  /* Optimal total points per iter = (dimension + 1) / 2
   * Optimal (per processor) points per iter = (dimension + 1) / 2 / size >= 1 <=> size <= (dimension + 1)/2 */
  int max_size = (dimension + 1)/2;
  size_used = std::min(size, max_size);
  
  this->rank = rank;
  int work_ranks[size_used];
  for(int i = 0; i < size_used; i++) {
    work_ranks[i] = i;
  }
  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group work_group;
  MPI_Group_incl(world_group, size_used, work_ranks, &work_group);
  MPI_Group_free(&world_group);
  MPI_Comm_create_group(MPI_COMM_WORLD, work_group, 0, &work_comm);
  MPI_Group_free(&work_group);
  
  global_best_par = new double[dimension];
  
  if(work_comm != MPI_COMM_NULL) {
    points_per_iter = max_size / size_used;
    
    /* Determine how many points are on the given processor, and their global
     * indices. Based off this index update with step size. */

    points_on_proc = (dimension + 1) / size_used;
    if ((dimension + 1) % size_used > rank)
      points_on_proc++;
    
    indices = new int[points_on_proc];
    for (int i = 0; i < points_on_proc; i++) {
      indices[i] = i;
    }
    
    simplex = new double[dimension * points_on_proc];
    create_simplex(&guess[0]);
    
    convergence = 0;
    this->obj_function = obj_function;
    this->obj_function_arg = obj_function_arg;
    M = new double[dimension];
    obj_function_results = new double[points_on_proc];
    AR = new double[dimension];
    AE = new double[dimension];
    AC = new double[dimension];
    updated = false;
    alpha = ALPHA;
    beta = BETA;
    gamma = GAMMA;
    tau = TAU;
    #if DEBUG
      logfile.open("log"+std::to_string(rank)+".txt");
    #endif
  }
}

DistParNelderMead::~DistParNelderMead() {
  delete global_best_par;
  
  if(work_comm != MPI_COMM_NULL) {
  	delete indices;
  	delete simplex;
  	delete M;
  	delete obj_function_results;
  	delete AR;
  	delete AE;
  	delete AC;
  	MPI_Comm_free(&work_comm);
  }
}

Rcpp::NumericVector DistParNelderMead::C_Array_to_R_Vec(double *vector, int dimension) {
  Rcpp::NumericVector R_Vec(vector, vector + dimension);
  return(R_Vec);
}

double DistParNelderMead::point_value(double* point) {
  Rcpp::NumericVector res = (*obj_function)(C_Array_to_R_Vec(point, dimension), obj_function_arg);
  return res[0];
}

Rcpp::List DistParNelderMead::solve(int max_iterations) {
  if(work_comm != MPI_COMM_NULL) {
  	// Compute objective function
  	for (int i = 0; i < points_on_proc; i++) {
  		obj_function_results[i] = point_value(&SIMPLEX(i, 0));
  	}
  	sort_simplex(); // Sort the simplex
    #if DEBUG
        print_simplex_with_value(); // Debugging purpose
    #endif
    global_best(work_comm);
  	
  	int iter = 0;
    double no_improve_tol = 10e-6;
    int no_improve_cnt = 0;
    int no_improve_max = 10;
  
    while(no_improve_cnt < no_improve_max && (max_iterations <= 0 || iter * size_used < max_iterations)) {
      #if DEBUG
  	    logfile << "Iter " << iter << "\n";
      #endif
  		current_point = points_on_proc - (iter % points_per_iter) - 1;
  
  		if (iter % points_per_iter == 0) {
  			centroid();
  		}
  		reflection(); //Compute reflection
  		fAR = point_value(AR); //Evaluate reflection

      #if DEBUG
    		logfile << "After reflection\n";
    		for(int i = 0; i < dimension; i++) {
    		  logfile << AR[i] << " ";
    		}
    		logfile << " ; Function value: " << fAR <<  "\n";
      #endif
  		
  		//Case 1
  		if (fAR < best) {
  			expansion();
  			fAE = point_value(AE); //Evaluate expansion
        #if DEBUG
    			logfile << "After expansion\n";
    			for(int i = 0; i < dimension; i++) {
    			  logfile << AE[i] << " ";
    			}
    			logfile << " ; Function value: " << fAE <<  "\n";
        #endif
    		//If expansion is better, use that
  			if (fAE < fAR) {
          #if DEBUG
  				  logfile << "expansion best: " << fAE << "\n";
          #endif
  				update(AE, current_point);
  				obj_function_results[indices[current_point]] = fAE;
  			} else { //otherwise use reflection
          #if DEBUG
  				  logfile << "reflection best\n";
          #endif
  				update(AR, current_point);
  				obj_function_results[indices[current_point]] = fAR;
  			}
  		//Case 2
  		} else if (fAR < obj_function_results[indices[current_point - 1]]) {
        #if DEBUG
  			  logfile << "AR better than next worst\n";
        #endif
  			update(AR, current_point);
  			obj_function_results[indices[current_point]] = fAR;
  		//Case 3
  		} else {
  			contraction();
  			fAC = point_value(AC); //Evaluate contraction
        #if DEBUG
    			logfile << "After contraction\n";
    			for(int i = 0; i < dimension;i++) {
    			  logfile << AC[i] << " ";
    			}
    			logfile << " ; Function value: " << fAC <<  "\n";
        #endif
  			// If Contraction is better, use it
  			if (fAC < std::min(obj_function_results[indices[current_point]], fAR)) {
          #if DEBUG
  				  logfile << "contraction better than worst\n";
          #endif
  				update(AC, current_point);
  				obj_function_results[indices[current_point]] = fAC;
  			} else { //Otherwise, minimize
  				if (fAR < obj_function_results[indices[current_point]]) {
            #if DEBUG
  					  logfile << "reflection better than worst\n";
            #endif
  					memmove(&SIMPLEX(current_point, 0), AR, dimension * sizeof(double));
  					obj_function_results[indices[current_point]] = fAR;
  				} //Else, not updated so don't bother
  			}
  		}
  		if ((iter % points_per_iter) == points_per_iter - 1) {
  			int global_updated = 0;
  			MPI_Allreduce(&updated, &global_updated, 1, MPI_INT, MPI_SUM,
                   work_comm);
  			if (!global_updated) { //not one processor had an update, minimize
  				minimize();
  				//Re-eval all of the points
  				for (int i = 0; i < points_on_proc; i++) {
  					obj_function_results[indices[i]] = point_value(&SIMPLEX(i, 0));
  				}
  			}
  			sort_simplex(); //Sort the simplex
        #if DEBUG
          	print_simplex_with_value();
        #endif
  			//Update global min on all processors
  			prev_best = best;
  			global_best(work_comm);
  			#if DEBUG
  			logfile << "global_best_value: " <<  best << "\n";
        #endif
  			if(best < prev_best - no_improve_tol) {
          no_improve_cnt = 0;
  			} else {
  			  no_improve_cnt += 1;
  			}
  			updated = 0;
  		}
  		iter++;
  	}
    if(rank == 0) {
      #if DEBUG
        Rcpp::Rcout << "Total Iterations: " << iter * size_used  << "\n";
      #endif
      if(max_iterations > 0 && iter * size_used >= max_iterations) {
        convergence = 1; // Iteration limit has been reached. Solution may not be correct.
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  global_best(MPI_COMM_WORLD);
  #if DEBUG
    logfile << "global_best_value: " <<  best << "\n";
    logfile.close();
  #endif
  Rcpp::NumericVector res_par = C_Array_to_R_Vec(global_best_par, dimension);
  MPI_Bcast(&convergence, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("par")=res_par, Rcpp::Named("value")=best, 
                                      Rcpp::Named("convergence")=convergence, Rcpp::Named("rank")=rank);
  return res;
}

void DistParNelderMead::create_simplex(double *guess) {
  /* e.g. 9 dimension, 10 (dimension + 1) points, 4 size_used
   * | Proc 0| Proc 1| Proc 2| Proc 3|
   * | 0 1 2 | 3 4 5 | 6 7   | 8 9   |
   * Let A be (dimension + 1) / size_used
   *     B be (dimension + 1) % size_used 
   * "More" processors (Proc 0 and 1) have (A+1) points, globalFirstIndex = (A+1)*rank = A*rank + rank when rank < B
   * "Less" processors (Proc 2 and 3) have A points, globalFirstIndex = (A+1)*B + A*(rank - B) = AB + B + A*rank - AB = A*rank + B when rank >= B */
  int globalFirstIndex = rank * ((dimension + 1) / size_used)
  + std::min((dimension + 1) % size_used, rank);
  
  /* Step size follows https://svn.r-project.org/R/trunk/src/appl/optim.c nmmin code,
   * which follows J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition */
  double step = 0;
  double step_cand;
  for(int i = 0; i < dimension; i++) {
    step_cand = 0.1*fabs(guess[i]);
    if(step_cand > step)
      step = step_cand;
  }
  if(step == 0) step = 0.1;
  
  for (int i = 0; i < points_on_proc; i++) {
    for (int j = 0; j < dimension; j++) {
      SIMPLEX(i, j) = guess[j];
      if (globalFirstIndex + i == j + 1)
        SIMPLEX(i, j) += step;
    }
  }
}


void DistParNelderMead::update(double *vector, int index) {
	if (!updated) { //only need to check if not already updated
		for (int i = 0; i < dimension; i++) {
			if (vector[i] != SIMPLEX(index, i)) {
				updated = 1;
				break;
			}
		}
	}
	if (updated) { //might be a new vector, copy it in
		memmove(&SIMPLEX(index, 0), vector, dimension * sizeof(double));
	}
}

void DistParNelderMead::centroid() {
	for (int i = 0; i < dimension; i++)
		M[i] = 0.0;
	for (int i = 0; i < (points_on_proc - points_per_iter); i++) {
		for (int j = 0; j < dimension; j++) {
			M[j] += SIMPLEX(i, j);
			//Divide after. Possible overflow for large obj function values!
		}
	}
	for (int i = 0; i < dimension; i++) {
		M[i] /= (dimension + 1 - points_per_iter); //Divide from earlier, then compute
	}
	// Reduce into AR, then swap pointers.
	MPI_Allreduce(M, AR, dimension, MPI_DOUBLE, MPI_SUM, work_comm);
	double *swap = M;
	M = AR;
	AR = swap;
}

void DistParNelderMead::reflection() {
	for (int i = 0; i < dimension; i++) {
		AR[i] = M[i] + alpha * (M[i] - SIMPLEX(current_point, i));
	}
}

void DistParNelderMead::expansion() {
	for (int i = 0; i < dimension; i++) {
		AE[i] = AR[i] + gamma * (AR[i] - M[i]);
	}
}

void DistParNelderMead::contraction() {
	double *ATilda;
	if (fAR < obj_function_results[indices[current_point]]) {
		ATilda = AR;
	} else {
		ATilda = &SIMPLEX(current_point, 0);
	}
	for (int i = 0; i < dimension; i++) {
		AC[i] = M[i] + beta * (ATilda[i] - M[i]);
	}
}

void DistParNelderMead::global_best(MPI_Comm comm) {
  struct {
    double val;
    int rank;
  } my_best_struct, global_best_struct;
  
  if(comm == work_comm || work_comm != MPI_COMM_NULL) {
    my_best_struct.val = obj_function_results[indices[0]];
  } else {
    my_best_struct.val = DBL_MAX;
  }
  my_best_struct.rank = rank;
  
  #if DEBUG
    if(work_comm != MPI_COMM_NULL) {
      logfile << "My best value: " << obj_function_results[indices[0]] << "\n";
      logfile << "My best par:\n";
      for(int i = 0; i < dimension ; i++) {
        logfile << SIMPLEX(0, i) << " ";
      }
      logfile << "obj_function_results:\n";
      for(int i = 0; i < points_on_proc ; i++) {
        logfile << obj_function_results[indices[i]] << " ";
      }
      logfile << "\n";
    }
  #endif
	
	MPI_Allreduce(&my_best_struct, &global_best_struct, 1, MPI_DOUBLE_INT, MPI_MINLOC,
               comm);
	if (rank == global_best_struct.rank) {
		memmove(global_best_par, &SIMPLEX(0, 0), dimension * sizeof(double));
	  best = obj_function_results[indices[0]];
	}
	MPI_Bcast(global_best_par, dimension, MPI_DOUBLE, global_best_struct.rank,
           comm);
	MPI_Bcast(&best, 1, MPI_DOUBLE, global_best_struct.rank,
           comm);
}

void DistParNelderMead::minimize() {
  global_best(work_comm);
	for (int i = 0; i < points_on_proc; i++) {
		daxpy(&SIMPLEX(i, 0), tau, global_best_par, (1.0 - tau), &SIMPLEX(i, 0),
				dimension);
	}
}

// result = scalar1*a + scalar2*b
void DistParNelderMead::daxpy(double *result, double scalar1, double *a,
		double scalar2, double *b, int length) {
	for (int i = 0; i < length; i++) {
		result[i] = scalar1 * a[i] + scalar2 * b[i];
	}
}

//Debugging purposes
#if DEBUG
  void DistParNelderMead::print_simplex() {
  	for (int i = 0; i < points_on_proc; i++) {
  		for (int j = 0; j < dimension; j++) {
  			Rcpp::Rcout << SIMPLEX(i, j) << " ";
  		}
  		Rcpp::Rcout << "\n";
  	}
  	Rcpp::Rcout << "\n";
  }
  void DistParNelderMead::print_simplex_with_value() {
    logfile << "From rank " << rank << "\n";
    logfile << "==================" << "\n";
    for (int i = 0; i < points_on_proc; i++) {
      for (int j = 0; j < dimension; j++) {
        logfile << SIMPLEX(i, j) << " ";
      }
      logfile << " ; Function value: " << point_value(&SIMPLEX(i, 0)) << "\n";
    }
    logfile << "\n";
  }
#endif

void DistParNelderMead::sort_simplex() {
	std::sort(indices, indices + points_on_proc,
			IndexSorter(obj_function_results));
}
