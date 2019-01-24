#!/usr/bin/Rscript
#  R/quadratic.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.24.2019
library(mvtnorm)
source('R/kern_funcs.R')


## Test that SOB
#' The quadratic function with linear and constant terms 0.
f <- function(x, A) drop(t(x) %*% A %*% x)

# Params
N <- 100
R <- 2
P <- 3

# Generate the target matrix
L <- matrix(rnorm(R*P), nrow = P)
A <- tcrossprod(L)

# Generate data
X <- matrix(runif(N * P), nrow = N)
y <- apply(X, 1, f, A = A)
# Center/scale it too
y <- (y - mean(y)) / sd(y)

init <- rnorm((R*P+P+P-R+1), 0, 1e-2)
init[(R*P+1):(R*P+P+P-R+1)] <- abs(init[(R*P+1):(R*P+P+P-R+1)])

lower_val <- 1e-4# How small are nonnegative quantities allowed to go?
ret <- optim(init, nll_wrapper, method = 'L-BFGS-B', 
      lower = c(rep(-Inf, R*P), rep(lower_val, (P+P-R+1))), y = y, X = X, R = R,
      control = list(trace = 5))
