# Run in terminal: mpirun -np 8 Rscript test.R

library(RDistParNelderMead)

compare_RDistParNelderMead_optim <- function(guess, fnc, list_args, correct_min, test_name) {
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
    start <- Sys.time()
    optim_res <- optim(guess, fnc, list_args=list_args, control=list(maxit=.Machine$integer.max))
    end <- Sys.time()
    print("optim result:")
    print(optim_res)
    print(end - start)
    if(abs(correct_min - optim_res$value)  < 1e-2) {
      print("Correct min")
    } else {
      print("Incorrect min")
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

compare_RDistParNelderMead_optim(rep(1, 100), objFunction1, list(0), 0, 1)
compare_RDistParNelderMead_optim(rep(1, 100), objFunction1, list(1), 0, 2)
compare_RDistParNelderMead_optim(rep(1, 100), objFunction2, list(0), 0, 3)
compare_RDistParNelderMead_optim(rep(1, 100), objFunction2, list(1), 0, 4)
RMPI_Finalize()
