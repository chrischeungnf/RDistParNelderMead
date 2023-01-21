# Run in terminal: mpirun -np 8 Rscript memory.R

library(RDistParNelderMead)
library(optimParallel)
setDefaultCluster(cl=makeCluster(8))

objFunction1 <- function(points, list_args) {
  shift <- list_args[[1]]
  return(sum((points + shift)^2/length(points)))
}

guess <- rep(1, 10000)
print("RDistParNelderMead:")
RDistParNelderMead_res <- RDistParNelderMead(guess, objFunction1, list(0))
RMPI_Finalize()
if(RDistParNelderMead_res$rank == 0) {
  print(RDistParNelderMead_res$value)
  print("optim using L-BFGS-B, serial:")
  optim_res <- optim(guess, objFunction1, list_args=list(0), method="L-BFGS-B", control=list(maxit=.Machine$integer.max))
  print(optim_res$value)
  print("optim using L-BFGS-B, parallel:")
  optim_parllel_res <- optimParallel(guess, objFunction1, list_args=list(0), control=list(maxit=.Machine$integer.max))
  print(optim_parllel_res$value)
}
