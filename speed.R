# Run in terminal: mpirun -np 8 Rscript speed.R

library(RDistParNelderMead)
library(tidyr)
library(ggplot2)

compare_speed <- function(guess, fnc, list_args, test_name) {
  power_vec <- seq(1, 4, 0.2)
  dim_vec <- as.integer(10^power_vec)
  RDistParNelderMead_duration <- optim_serial_duration <- numeric(length(dim_vec))
  
  for(i in seq_along(dim_vec)) {
    guess_vec <- rep(guess, dim_vec[i])
    start <- Sys.time()
    RDistParNelderMead_res <- RDistParNelderMead(guess_vec, fnc, list_args)
    end <- Sys.time()
    if(RDistParNelderMead_res$rank == 0) {
      RDistParNelderMead_duration[i] <- as.numeric(end - start)
      
      start <- Sys.time()
      optim(guess_vec, fnc, list_args=list_args, method="L-BFGS-B", control=list(maxit=.Machine$integer.max))
      end <- Sys.time()
      optim_serial_duration[i] <- as.numeric(end - start)
    }
  }
  if(RDistParNelderMead_res$rank == 0) {
    speed_df <- data.frame(dimension=dim_vec, RDistParNelderMead=RDistParNelderMead_duration,
                           optim_serial=optim_serial_duration)
    speed_df_long <- speed_df %>% pivot_longer(!dimension, names_to = "method", values_to = "seconds")
    speed_plot <- ggplot(data = speed_df_long, aes(x=dimension, y=seconds)) + geom_line(aes(colour=method))
    ggsave(paste0(test_name, ".png"), speed_plot)
    
    speed_df$RDistParNelderMead_best <- speed_df$RDistParNelderMead < speed_df$optim_serial
    print(paste("Test", test_name, "results:"))
    print("==============")
    print(speed_df)
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

compare_speed(1, objFunction1, list(0), "S1")
compare_speed(1, objFunction2, list(0), "S2")
RMPI_Finalize()
