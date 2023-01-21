objFunction1 <- function(points, list_args) {
  shift <- list_args[[1]]
  return(sum((points + shift)^2/length(points)))
}

library(RDistParNelderMead)
shift <- 1
res <- RDistParNelderMead(rep(1, 100), objFunction1, list(shift))
if(res$rank == 0) {
    print(res)
}
RMPI_Finalize()

