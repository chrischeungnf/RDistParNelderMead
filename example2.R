objFunction1_nopara <- function(points, list_args) {
  return(sum(points^2/length(points)))
}

library(RDistParNelderMead)
res <- RDistParNelderMead(rep(1, 100), objFunction1_nopara)
if(res$rank == 0) {
    print(res)
}
RMPI_Finalize()

