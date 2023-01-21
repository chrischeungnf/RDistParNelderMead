# RDistParNelderMead
MPI-Based Distributed Memory Parallel Nelder-Mead Method for R

## Description
This project brings MPI-Based Distributed Memory Parallel Nelder-Mead Method proposed in Klein and Neira (2014) to R users. The problem this project solves is unconstrained minimization for high-dimensional functions.

## Installation
The current supported operating system is Linux. The Open MPI library is required for the code. If your operating system is Debian or its derivatives, you can install Open MPI library by:
```
sudo apt install libopenmpi-dev
```
Then, you can install devtools library on R if you have not already done so. On Debian or its derivatives, you can run:
```
sudo apt install r-cran-devtools
```
Finally, you can install RDistParNelderMead library on R using devtools:
```
library(devtools)
install_github("chrischeungnf/RDistParNelderMead")
```

## Usage
Wrap your R code which uses RDistParNelderMead library in a file. Then run the file as follows:
```
mpirun -np x Rscript y.R
```
where x is the number of processors you want to use and y.R is your R file.

## Example
### Example 1 (example1.R)
For example, if you want to minimize the following objective function:
$$f_1(x, shift) = \sum_{i=1}^n \frac{(x_i + shift)^2}{n}$$

You can first set up the objective function in R:
```
objFunction1 <- function(points, list_args) {
  shift <- list_args[[1]]
  return(sum((points + shift)^2/length(points)))
}
```

Note that RDistParNelderMead assumes function passed to it have 2 arguments. The first argument accepts a numeric vector while second argument accepts a list of parameters for the function.

Then, you can minimize the function with $shift = 1$, $n = 100$ and initial guess $(1, 1, 1, ..., 1)$ by:
```
library(RDistParNelderMead)
shift <- 1
res <- RDistParNelderMead(rep(1, 100), objFunction1, list(shift))
if(res$rank == 0) {
    print(res)
}
```

After you have used RDistParNelderMead, terminate the MPI execution environment by: 
```
RMPI_Finalize()
```

After that, you can wrap the code in a R script and run it with x processors:
```
mpirun -np x Rscript example1.R
```

Result:
```
$par
  [1] -0.9995026 -0.9995026 -0.9990283 -0.9990283 -0.9990283 -0.9990283
  [7] -1.0046810 -0.9886285 -0.9948785 -1.0003472 -1.0007379 -1.0009332
 [13] -0.9966033 -0.9995026 -0.9995026 -0.9990283 -0.9990283 -0.9990283
 [19] -0.9990283 -1.0047542 -1.0011285 -1.0011285 -1.0011285 -1.0011285
 [25] -1.0011285 -0.9964080 -0.9995026 -0.9995026 -0.9990283 -0.9990283
 [31] -0.9990283 -0.9990283 -1.0047786 -1.0011285 -1.0011285 -1.0011285
 [37] -1.0011285 -1.0011285 -0.9962188 -0.9995026 -0.9995026 -0.9990283
 [43] -0.9990283 -0.9990283 -0.9990283 -1.0047786 -1.0011285 -1.0011285
 [49] -1.0011285 -1.0011285 -1.0011285 -0.9962158 -0.9995026 -0.9995026
 [55] -0.9990283 -0.9990283 -0.9990283 -0.9990283 -1.0047786 -1.0011285
 [61] -1.0011285 -1.0011285 -1.0011285 -1.0011285 -0.9990283 -0.9990283
 [67] -0.9990283 -0.9990283 -0.9990283 -0.9990283 -1.0011285 -1.0011285
 [73] -1.0011285 -1.0011285 -1.0011285 -1.0011285 -0.9990283 -0.9990283
 [79] -0.9990283 -0.9990283 -0.9990283 -0.9990283 -1.0011285 -1.0011285
 [85] -1.0011285 -1.0011285 -1.0011285 -1.0011285 -0.9990283 -0.9990283
 [91] -0.9990283 -0.9990283 -0.9990283 -0.9990283 -1.0011285 -1.0011285
 [97] -1.0011285 -1.0011285 -1.0011285 -1.0011285

$value
[1] 4.099012e-06

$convergence
[1] 0

$rank
[1] 0
```

RDistParNelderMead returns a list containing par (optimal point), value (minimized value), convergence (0 for convergence; 1 for iteration limit has been reached, which indicates optimal point may not be correct) and rank (Processor ID).

### Example 2 (example2.R)
If your objective function does not have parameter like:
$$f_1(x) = \sum_{i=1}^n \frac{(x_i)^2}{n}$$

Your R function still needs to accept two arguments:
```
objFunction1_nopara <- function(points, list_args) {
  return(sum(points^2/length(points)))
}
```
But, you can minimize the function by:
```
library(RDistParNelderMead)
res <- RDistParNelderMead(rep(1, 100), objFunction1_nopara)
```

Result:
```
$par
  [1] -2.784329e-05 -7.322806e-05 -1.864472e-05 -7.226434e-05 -1.437036e-05
  [6] -2.584871e-05  1.343947e-05  9.465804e-05  6.354349e-05 -3.691190e-05
 [11]  1.266141e-04 -1.644155e-05 -6.569800e-05 -2.784329e-05 -7.322806e-05
 [16] -1.864472e-05 -7.226434e-05 -1.437036e-05 -2.584871e-05  1.343947e-05
 [21]  9.465804e-05  6.354349e-05 -3.691190e-05  1.266141e-04 -1.644155e-05
 [26] -6.569800e-05 -2.784329e-05 -7.322806e-05 -1.864472e-05 -7.226434e-05
 [31] -1.437036e-05 -2.584871e-05  1.343947e-05  9.465804e-05  6.354349e-05
 [36] -3.691190e-05  1.266141e-04 -1.644155e-05 -6.569800e-05 -2.784329e-05
 [41] -7.322806e-05 -1.864472e-05 -7.226434e-05 -1.437036e-05 -2.584871e-05
 [46]  1.343947e-05  9.465804e-05  6.354349e-05 -3.691190e-05  1.266141e-04
 [51] -1.644155e-05 -6.569800e-05 -2.784329e-05 -7.322806e-05 -1.864472e-05
 [56] -7.226434e-05 -1.437036e-05 -2.584871e-05  1.343947e-05  9.465804e-05
 [61]  6.354349e-05 -3.691190e-05  1.266141e-04 -1.644155e-05 -5.054524e-05
 [66] -4.316041e-05 -4.316041e-05  3.416776e-04 -4.894735e-05 -5.365932e-05
 [71] -7.415574e-07  3.050697e-06 -6.836273e-07  4.204247e-05 -3.425323e-06
 [76]  5.582431e-05 -5.054524e-05 -4.316041e-05 -4.316041e-05 -4.894735e-05
 [81] -4.894735e-05 -5.365932e-05 -7.415574e-07  3.050697e-06 -6.836273e-07
 [86]  4.204247e-05 -3.425323e-06  5.582431e-05 -5.054524e-05 -4.316041e-05
 [91] -4.316041e-05 -4.894735e-05 -4.894735e-05 -5.365932e-05 -7.415574e-07
 [96]  3.050697e-06 -6.836273e-07  4.204247e-05 -3.425323e-06  5.582431e-05

$value
[1] 4.053075e-09

$convergence
[1] 0

$rank
[1] 0
```

### Other Examples
For other examples, you can refer to test.R on this repo.

## Citation
Klein, K., and J. Neira (2014). Nelder-Mead Simplex Optimization Routine for Large-Scale Problems: A Distributed Memory Implementation. *Computational Economics, 43*(4), 447–461. <https://doi.org/10.1007/s10614-013-9377-8>
