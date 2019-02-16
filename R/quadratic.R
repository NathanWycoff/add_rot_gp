#!/usr/bin/Rscript
#  R/quadratic.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.24.2019
library(mvtnorm)
library(activegp)
library(dr)
library(hetGP)
library(DiceDesign)
library(Rcpp)
library(foreach)
library(doParallel)
source('R/kern_funcs.R')
sourceCpp('src/kern_funcs.cpp')#, rebuild = T, verbose = T)

## Test that SOB
set.seed(1234)
#' The quadratic function with linear and constant terms 0.
f <- function(x, A) drop(t(x) %*% A %*% x)

# Params
n_cores <- 20
registerDoParallel(n_cores)
sets <- 3
R <- 1
van_P <- 10 # If dimension is at least this much, skip vanilla.
Ps <- c(2, 10, 50)
Ns <- c(30, 100, 300)
NN <- 1000# num of prediction locations
iterss <- c(40, 40, 40)# How many sims to run?

results <- list()
save_file <- paste(file.path('data', 'sims', 'quadratic'), as.numeric(Sys.time()), '.RData', sep = '')


lhs_seeds <- lapply(1:sets, function(set) sample(1e6,iterss[set]))

for (set in 1:sets) {
    cat("Set:", set, '-------------------\n')

    # Get settings
    P <- Ps[set]
    N <- Ns[set]
    iters <- iterss[set]

    # Storage for results.
    if (P < van_P) {
        nm_1 <- 10
        efforts <- matrix(NA, nrow = iters, ncol = nm_1)
        nlls <- matrix(NA, nrow = iters, ncol = nm_1)
        is_errs <- matrix(NA, nrow = iters, ncol = nm_1)
        oos_errs <- matrix(NA, nrow = iters, ncol = nm_1)
        sub_errs <- matrix(NA, nrow = iters, ncol = nm_1)

        cnames_1 <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S", "VAN", "VAN_S", "ADD", "ADD_S")
        eval_subspace <- c(T,T,T,T,T,T,F,F,F,F)#Does a method try to estimate an important subspace?
        colnames(efforts) <- cnames_1
        colnames(nlls) <- cnames_1
        colnames(is_errs) <- cnames_1
        colnames(oos_errs) <- cnames_1
        colnames(sub_errs) <- cnames_1
    } else {
        nm_2 <- 8
        efforts <- matrix(NA, nrow = iters, ncol = nm_2)
        nlls <- matrix(NA, nrow = iters, ncol = nm_2)
        is_errs <- matrix(NA, nrow = iters, ncol = nm_2)
        oos_errs <- matrix(NA, nrow = iters, ncol = nm_2)
        sub_errs <- matrix(NA, nrow = iters, ncol = nm_2)

        cnames_2 <- c("MF", "MF_S", "MF_P", "MF_SP", "LR", "LR_S", "ADD", "ADD_S")
        eval_subspace <- c(T,T,T,T,T,T,F,F)
        colnames(efforts) <- cnames_2
        colnames(nlls) <- cnames_2
        colnames(is_errs) <- cnames_2
        colnames(oos_errs) <- cnames_2
        colnames(sub_errs) <- cnames_2
    }

    null_mse <- rep(NA, iters)
    hetgp_mse <- rep(NA, iters)
    true_var <- rep(NA, iters)


    #for (iter in 1:iters) {
    parret <- foreach(iter = 1:iters) %dopar% {
        cat("Iter:", iter, '\n')
        cat("Seed:", lhs_seeds[[sets]][iter], '\n')
        set.seed(lhs_seeds[[sets]][iter])
        # Generate the target matrix
        B <- qr.Q(qr(matrix(rnorm(P*P), nrow = P)))
        LAMBDA <- diag(10^(1:(-P+2)))
        A <- B %*% LAMBDA %*% t(B)
        L <- B[,1:R, drop = FALSE]
        #L <- matrix(rnorm(R*P), nrow = P)
        #A <- tcrossprod(L) / 10

        # Generate data
        #X <- matrix(runif(N * P), nrow = N)
        X <- lhsDesign(N,P, seed = lhs_seeds[[sets]][iter])$design

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
        ret_ad  <- gp_mf(y, X, R = 0, has_tail = TRUE)
        ret_ad_s  <- gp_mf(y, X, R = 0, has_tail = TRUE, smart_init = FALSE)
        if (P < 10) {
            ret_van  <- gp_mf(y, X, R = P)
            ret_van_s  <- gp_mf(y, X, R = P, smart_init = FALSE)
            rets <- list(ret_mf, ret_mf_s, ret_mf_p, ret_mf_sp, 
                         ret_lr, ret_lr_s, ret_van, ret_van_s, ret_ad, ret_ad_s)
        } else {
            rets <- list(ret_mf, ret_mf_s, ret_mf_p, ret_mf_sp, 
                         ret_lr, ret_lr_s, ret_ad, ret_ad_s)
        }

        # hetGP model
        ret_gp <- mleHomGP(X, y, lower = rep(1e-2, P), upper = rep(P,P),
                        noiseControl = list(g_bounds = c(1e-6, 1e-6)), covtype = 'Gaussian')

        # Evaluate 
        res <- list(sapply(rets, function(mod) mod$tt[3]),
          sapply(rets, function(mod) mod$optim_ret$val),
          sapply(rets, function(mod) mean((predict(mod, X) - y)^2)),
          sapply(rets, function(mod) mean((predict(mod, XX) - yy)^2)),
          sapply(1:length(rets), function(m) ifelse(eval_subspace[m], subspace_dist(rets[[m]]$U[,1:R],L), NA)))

        # Baseline metrics: Standard GP and null model
        res <- c(res, mean((mean(y) - yy)^2),
                 mean((predict(ret_gp, XX)$mean - yy)^2),
                 var(yy))

        return(res)
    }

    for (iter in 1:iters) {
        efforts[iter,] <- parret[[iter]][[1]]
        nlls[iter,] <- parret[[iter]][[2]]
        is_errs[iter,] <- parret[[iter]][[3]]
        oos_errs[iter,] <- parret[[iter]][[4]]
        sub_errs[iter,] <- parret[[iter]][[5]]

        null_mse[iter] <- parret[[iter]][[6]]
        hetgp_mse[iter] <- parret[[iter]][[7]]
        true_var[iter] <- parret[[iter]][[8]]
    }

    results[[set]] <- list(efforts, nlls, is_errs, oos_errs, sub_errs, null_mse, hetgp_mse, true_var)

    save(results, R, Ps, Ns, NN, iterss,
         file = save_file)
}
stopImplicitCluster()
