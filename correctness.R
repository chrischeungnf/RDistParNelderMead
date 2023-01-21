# Run in terminal: mpirun -np 8 Rscript correctness.R

library(RDistParNelderMead)
library(optimParallel)
setDefaultCluster(cl=makeCluster(8))

compare_RDistParNelderMead_optim <- function(guess, fnc, list_args, optim_methods, correct_min, test_name) {
  start <- Sys.time()
  RDistParNelderMead_res <- RDistParNelderMead(guess, fnc, list_args)
  end <- Sys.time()
  if(RDistParNelderMead_res$rank == 0) {
    print(paste("Test", test_name, "results:"))
    print("==============")
    print("RDistParNelderMead result:")
    print(RDistParNelderMead_res)
    print(end - start)
    if(abs(correct_min - RDistParNelderMead_res$value)  < 1e-2) {
      print("Correct min")
    } else {
      print("Incorrect min")
    }
    print("--------------")
    for(method in optim_methods) {
      start <- Sys.time()
      if(method == "parallel-L-BFGS-B") {
        optim_res <- optimParallel(guess, fnc, list_args=list_args, control=list(maxit=.Machine$integer.max))
      } else {
        optim_res <- optim(guess, fnc, list_args=list_args, method=method, control=list(maxit=.Machine$integer.max))
      }
      end <- Sys.time()
      print(paste(method, "result:"))
      print(optim_res)
      print(end - start)
      if(abs(correct_min - optim_res$value)  < 1e-2) {
        print("Correct min")
      } else {
        print("Incorrect min")
      }
      print("--------------")
    }
  }
}

objFunction1 <- function(points, list_args) {
  shift <- list_args[[1]]
  return(sum((points + shift)^2/length(points)))
}

objFunction2 <- function(points, list_args) {
  shift <- list_args[[1]]
  return(sum(abs(points + shift)/length(points)))
}

compare_RDistParNelderMead_optim(rep(1, 100), objFunction1, list(0), c("Nelder-Mead", "L-BFGS-B", "parallel-L-BFGS-B"), 0, "C1")
compare_RDistParNelderMead_optim(rep(1, 100), objFunction2, list(0), c("Nelder-Mead", "L-BFGS-B", "parallel-L-BFGS-B"), 0, "C2")
compare_RDistParNelderMead_optim(rep(1, 100), objFunction1, list(1), c("Nelder-Mead", "L-BFGS-B", "parallel-L-BFGS-B"), 0, "C3")
compare_RDistParNelderMead_optim(rep(1, 100), objFunction2, list(1), c("Nelder-Mead", "L-BFGS-B", "parallel-L-BFGS-B"), 0, "C4")
RMPI_Finalize()
