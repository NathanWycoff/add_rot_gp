#!/usr/bin/Rscript
#  R/quadratic.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.24.2019
library(mvtnorm)
library(activegp)
library(dr)
library(hetGP)
library(DiceDesign)
library(Rcpp)
source('R/kern_funcs.R')
sourceCpp('src/kern_funcs.cpp')

## Test that SOB
set.seed(1234)
#' The quadratic function with linear and constant terms 0.
f <- function(x, A) drop(t(x) %*% A %*% x)

# Params
sets <- 3
R <- 1
van_P <- 10 # If dimension is at least this much, skip vanilla.
Ps <- c(2, 10, 50)
Ns <- c(30, 100, 300)
NN <- 1000# num of prediction locations
iterss <- c(30, 30, 5)# How many sims to run?

results <- list()
save_file <- paste(file.path('data', 'sims', 'quadratic'), as.numeric(Sys.time()), '.RData', sep = '')

for (set in 1:sets) {
    cat("Set:", set, '\n')

    # Get settings
    P <- Ps[set]
    N <- Ns[set]
    iters <- iterss[set]

    # Storage for results.
    if (P < van_P) {
        efforts <- matrix(NA, nrow = iters, ncol = 8)
        nlls <- matrix(NA, nrow = iters, ncol = 8)
        is_errs <- matrix(NA, nrow = iters, ncol = 8)
        oos_errs <- matrix(NA, nrow = iters, ncol = 8)
        sub_errs <- matrix(NA, nrow = iters, ncol = 8)

        colnames(efforts) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S", "VAN", "VAN_S")
        colnames(nlls) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S", "VAN", "VAN_S")
        colnames(is_errs) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S", "VAN", "VAN_S")
        colnames(oos_errs) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S", "VAN", "VAN_S")
        colnames(sub_errs) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S", "VAN", "VAN_S")
    } else {
        efforts <- matrix(NA, nrow = iters, ncol = 6)
        nlls <- matrix(NA, nrow = iters, ncol = 6)
        is_errs <- matrix(NA, nrow = iters, ncol = 6)
        oos_errs <- matrix(NA, nrow = iters, ncol = 6)
        sub_errs <- matrix(NA, nrow = iters, ncol = 6)

        colnames(efforts) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S")
        colnames(nlls) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S")
        colnames(is_errs) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S")
        colnames(oos_errs) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S")
        colnames(sub_errs) <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S")
    }

    null_mse <- rep(NA, iters)
    hetgp_mse <- rep(NA, iters)
    true_var <- rep(NA, iters)


    for (iter in 1:iters) {
        cat("Iter:", iter, '\n')
        # Generate the target matrix
        B <- qr.Q(qr(matrix(rnorm(P*P), nrow = P)))
        LAMBDA <- diag(10^(1:(-P+2)))
        A <- B %*% LAMBDA %*% t(B)
        L <- B[,1:R, drop = FALSE]
        #L <- matrix(rnorm(R*P), nrow = P)
        #A <- tcrossprod(L) / 10

        # Generate data
        #X <- matrix(runif(N * P), nrow = N)
        X <- lhsDesign(N,P)$design

        XX <- matrix(runif(NN * P, min = 0.1, max = 0.9), nrow = NN)
        y <- apply(X, 1, f, A = A)
        yy <- apply(XX, 1, f, A = A)
        # Center/scale it too, according to training set mean/var.
        y <- (y - mean(y)) / sd(y)
        yy <- (yy - mean(y)) / sd(y)

        # Fit our custom methods
        ret_mf  <- gp_mf(y, X, R = R, smart_init = FALSE)
        ret_mf_s  <- gp_mf(y, X, R = R)
        ret_mf_p  <- gp_mf(y, X, R = R, smart_init = FALSE, preinit = TRUE)
        ret_mf_sp  <- gp_mf(y, X, R = R, smart_init = TRUE, preinit = TRUE)
        ret_lr  <- gp_mf(y, X, R = R, has_tail = FALSE)
        ret_lr_s  <- gp_mf(y, X, R = R, has_tail = FALSE, smart_init = FALSE)
        if (P < 10) {
            ret_van  <- gp_mf(y, X, R = P)
            ret_van_s  <- gp_mf(y, X, R = P, smart_init = FALSE)
            rets <- list(ret_mf, ret_mf_s, ret_mf_p, ret_mf_sp, 
                         ret_lr, ret_lr_s, ret_van, ret_van_s)
        } else {
            rets <- list(ret_mf, ret_mf_s, ret_mf_p, ret_mf_sp, 
                         ret_lr, ret_lr_s)
        }

        # hetGP model
        ret_gp <- mleHomGP(X, y, lower = rep(1e-4, P), upper = rep(P,P),
                        noiseControl = list(g_bounds = c(1e-6, 1e-6)), covtype = 'Gaussian')

        # Evaluate 
        efforts[iter,] <- sapply(rets, function(mod) mod$tt[3])
        nlls[iter,] <- sapply(rets, function(mod) mod$optim_ret$val)
        is_errs[iter,] <- sapply(rets, function(mod) mean((predict(mod, X) - y)^2))
        oos_errs[iter,] <- sapply(rets, function(mod) mean((predict(mod, XX) - yy)^2))
        sub_errs[iter,] <- sapply(rets, function(mod) subspace_dist(mod$U[,1:R],L))

        # Baseline metrics: Standard GP and null model
        null_mse[iter] <- mean((mean(y) - yy)^2)
        hetgp_mse[iter] <- mean((predict(ret_gp, XX)$mean - yy)^2)
        true_var[iter] <- var(yy)
    }

    results[[set]] <- list(efforts, nlls, is_errs, oos_errs, sub_errs, null_mse, hetgp_mse, true_var)

    save(results, R, Ps, Ns, NN, iterss,
         file = save_file)
}
