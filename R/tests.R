#!/usr/bin/Rscript
#  R/tests.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.24.2019

library(microbenchmark)
library(mvtnorm)
library(activegp)
library(dr)
library(hetGP)
library(DiceDesign)
library(Rcpp)
source('R/kern_funcs.R')
sourceCpp('src/kern_funcs.cpp')

## Verify that Cpp agrees with R code.
# Kernel itself
P <- 4
R <- 2

rot <- matrix(rnorm(P*P), nrow = P)
U <- rot[,1:R]
if (R < P) {
    Vs <- lapply((R+1):P, function(p) rot[,p,drop=FALSE])
} else {
    Vs <- NULL
}

d <- rnorm(P, 0, 1e-2)
l <- abs(rnorm(P, 0, 1e-2))
gamma <- abs(rnorm(P-R+1,0,1e-2))

fast_k(d, gamma, R, TRUE)
fast_k_cpp(d, gamma, R, TRUE)
fast_k(d, gamma, R, FALSE)
fast_k_cpp(d, gamma, R, FALSE)

microbenchmark(fast_k(d, gamma, R, TRUE))
microbenchmark(fast_k_cpp(d, gamma, R, TRUE))

# Forming the covariance matrix
#' The quadratic function with linear and constant terms 0.
f <- function(x, A) drop(t(x) %*% A %*% x)

# Params
N <- 100
NN <- 10# num of prediction locations
R <- 2
P <- 6

# Generate the target matrix
L <- matrix(rnorm(R*P), nrow = P)
A <- tcrossprod(L) / 10

# Generate data
#X <- matrix(runif(N * P), nrow = N)
X <- lhsDesign(N,P)$design

XX <- matrix(runif(NN * P, min = 0.1, max = 0.9), nrow = NN)
y <- apply(X, 1, f, A = A)
yy <- apply(XX, 1, f, A = A)
# Center/scale it too, according to training set mean/var.
y <- (y - mean(y)) / sd(y)
yy <- (yy - mean(y)) / sd(y)

# Build the covariance matrix
rot <- matrix(rnorm(P*P), nrow = P)
U <- rot[,1:R]
if (R < P) {
    Vs <- lapply((R+1):P, function(p) rot[,p,drop=FALSE])
} else {
    Vs <- NULL
}

d <- rnorm(P, 0, 1e-2)
l <- abs(rnorm(P, 0, 1e-2))
gamma <- abs(rnorm(P-R+1,0,1e-2))

Xr <- X %*% rot
XXr <- XX %*% rot
    
form_cov(Xr, k = function(d) fast_k(d, gamma = gamma, R = R, has_tail = TRUE), g = 1e-6) - form_cov_cpp(Xr, gamma, R, TRUE, 1e-6)
form_cov(Xr, k = function(d) fast_k(d, gamma = gamma, R = R, has_tail = FALSE), g = 1e-6) - form_cov_cpp(Xr, gamma, R, FALSE, 1e-6)

form_cov(Xr, XXr, k = function(d) fast_k(d, gamma = gamma, R = R, has_tail = TRUE), g = 1e-6) - form_crosscov_cpp(Xr, XXr, gamma, R, TRUE, 1e-6)
form_cov(Xr, XXr, k = function(d) fast_k(d, gamma = gamma, R = R, has_tail = FALSE), g = 1e-6) - form_crosscov_cpp(Xr, XXr, gamma, R, FALSE, 1e-6)

## Test at NLL level.
microbenchmark(nll(y, X, R, U, Vs, l, gamma, g = 1e-6, has_tail = TRUE))

## Speedups, verify that speedups are still accurate
P <- 4
R <- 2

rot <- matrix(rnorm(P*P), nrow = P)
U <- rot[,1:R]
if (R < P) {
    Vs <- lapply((R+1):P, function(p) rot[,p,drop=FALSE])
} else {
    Vs <- NULL
}

d <- rnorm(P, 0, 1e-2)
l <- abs(rnorm(P, 0, 1e-2))
gamma <- abs(rnorm(P-R+1,0,1e-2))

dr <- rot %*% d
drs <- dr / sqrt(l)

microbenchmark(old_k(d, U, Vs, l, gamma))
microbenchmark(k(d, U, Vs, l, gamma))
microbenchmark(fast_k(d, gamma))
