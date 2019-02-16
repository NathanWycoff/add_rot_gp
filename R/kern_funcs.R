#!/usr/bin/Rscript
#  R/kern_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.24.2019

#' The kernel to be studied
#' @param d The difference between two points.
#' @param U The high fidelity subspace, or, to be precise, a matrix the columns of which contain an orthonormal basis therefor.
#' @param Vs A list containing unit vectors filling the space after U.
#' @param l A numeric vector of nonnegative lengthscales. Gives the variance if we were to consider this a Gaussian.
#' @param gamma The covariance or scaling factors for each component of the sum; should be of length P - R + 1. 
#' @param has_tail Do we include the additive tail?
k <- function(d, U, Vs, l, gamma, has_tail = TRUE) {
    if (is.null(dim(U))) {
        R <- 0
    } else {
        R <- ncol(U)
    }
    Pvs <- lapply(Vs, tcrossprod)
    corr <- 0
    #corr <- gamma[1] * exp(-0.5 * t(d) %*% U %*% diag(1/l[1:R]) %*% t(U) %*% d)
    if (R > 0) {
        corr <- corr + gamma[1] * exp(-0.5 * sum((t(d) %*% U %*% diag(1/sqrt(l[1:R]), nrow = R))^2))
    }
    if (R < P) {
        if (has_tail) {
            for (p in (R+1):P) {
                #corr <- corr + gamma[p-R+1] * exp(-0.5 * t(d) %*% Pvs[[p-R]] %*% d / l[p])
                #print(p)
                #print(p-R)
                #print(Vs)
                corr <- corr + gamma[p-R+1] * exp(-0.5 * sum((t(d) %*% Vs[[p-R]])^2) / l[p])
            }
        }
    }
    return(corr)
}

#' The kernel to be studied
#' @param d The difference between two points, it is expected to have already been rotated and scaled.
#' @param gamma The covariance or scaling factors for each component of the sum; should be of length P - R + 1
fast_k <- function(d, gamma, R, has_tail = TRUE) {
    P <- length(d)
    corr <- 0
    if (R > 0) {
    corr <- corr + gamma[1] * exp(-0.5 * sum(d[1:R]^2))
    }
    if (has_tail && R < P) {
        for (p in (R+1):P) {
            corr <- corr + gamma[p-R+1] * exp(-0.5 * (d[p]^2))
        }
    }
    return(corr)
}

#' Form GP Covariance structure. Nugget added to diagonal if X2 not given.
form_cov <- function(X1, X2 = NULL, k, g = 1e-6) {
    if (is.null(X2)) {
        K <- matrix(NA, nrow = nrow(X1), ncol = nrow(X1))
        for (n1 in 1:nrow(X1)) {
            for (n2 in n1:nrow(X1)) {
                K[n1,n2] <- K[n2,n1] <- k(X1[n1,] - X1[n2,])
            }
        }
        diag(K) <- diag(K) + g
    } else {
        K <- matrix(NA, nrow = nrow(X1), ncol = nrow(X2))
        for (n1 in 1:nrow(X1)) {
            for (n2 in 1:nrow(X2)) {
                K[n1,n2] <- k(X1[n1,] - X2[n2,])
            }
        }
    }
    return(K)
}

#' Negative GP Log Likelihood.
#' @param y The response, a numeric.vector.
#' @param X The design, a matrix with as many rows as observations.
#' @param U The high fidelity subspace, or, to be precise, a matrix the columns of which contain an orthonormal basis therefor.
#' @param Vs A list containing unit vectors filling the space after U.
#' @param l A numeric vector of nonnegative lengthscales. Gives the variance if we were to consider this a Gaussian.
#' @param gamma The covariance or scaling factors for each component of the sum; should be of length P - R + 1
#' @param g The error variance.
#' @param has_tail Do we include the additive tail?
nll <- function(y, X, R, U, Vs, l, gamma, g = 1e-6, has_tail = TRUE) {
    N <- nrow(X)
    P <- ncol(X)

    # Build the covariance matrix
    Xr <- X %*% do.call(cbind, c(list(U), Vs))
    Xrs <- Xr %*% diag(1/sqrt(l))
    #K <- form_cov(Xrs, k = function(d) fast_k(d, gamma = gamma, R = R, has_tail = has_tail), g = g)
    K <- form_cov_cpp(Xrs, gamma, R, has_tail, g)

    # Old func:
    #K <- form_cov(X, k = function(d) k(d, U, Vs, l, gamma, has_tail), g = g)
    #K <- matrix(NA, nrow = N, ncol = N)
    #for (n1 in 1:N) {
    #    for (n2 in n1:N) {
    #        K[n1,n2] <- K[n2, n1] <- k(X[n1,] - X[n2,], U, Vs, l, gamma, has_tail)
    #    }
    #}
    #diag(K) <- diag(K) + g

    # Evaluate loglik
    return(-dmvnorm(y, rep(0, N),  K, log = T))
}

#TODO: Use this in the nll_wrapper.
par_2_param <- function(par, P, R) {
    # Form parameters
    L <- matrix(par[1:(R*P)], nrow = P)
    l <- par[(R*P+1):(R*P + P)]
    gamma <- par[(R*P+P+1):(R*P+P+P-R+1)]

    # Build orthonormal basis.
    B <- qr.Q(qr(cbind(L, matrix(rnorm((P-R)*P), nrow = P))))
    if (R > 0) {
        U <- B[,1:R, drop = FALSE]
    } else {
        U <- numeric(0)
    }
    if (R < P) {
        Vs <- lapply((R+1):P, function(p) B[,p, drop = FALSE])
    } else {
        Vs <- NULL
    }

    return(list(U = U, Vs = Vs, l = l, gamma = gamma))
}

param_2_par <- function(U, l, gamma) {
    return(c(as.numeric(U), l, gamma))
}

#' A wrapper of the negative loglikelihood for use by optim
#' @param par is a vector, the first R*P elements give a matrix the range of which is a basis for the HF space. The next P parameters give lengthscales, and the P-R+1 parameters after that give covariances. As such, this vector should be of length (R*P+P+P-R+1). The first R*P are unconstrained, while the remaining should be bounded below by 0, or, preferably, a small constant.
nll_wrapper <- function(par, y, X, R, has_tail) {
    N <- nrow(X)
    P <- ncol(X)

    # Form parameters
    L <- matrix(par[1:(R*P)], nrow = P)
    l <- par[(R*P+1):(R*P + P)]
    gamma <- par[(R*P+P+1):(R*P+P+P-R+1)]

    # Build orthonormal basis.
    #TODO: Better init
    set.seed(123)
    B <- qr.Q(qr(cbind(L, matrix(rnorm((P-R)*P), nrow = P))))
    if (R > 0) {
        U <- B[,1:R, drop = FALSE]
    } else {
        U <- numeric(0)
    }
    if (R < P) {
        Vs <- lapply((R+1):P, function(p) B[,p, drop = FALSE])
    } else {
        Vs <- NULL
    }

    #print("U")
    #print(U)
    #print("l")
    #print(l)
    #print("gamma")
    #print(gamma)

    ret <- nll(y, X, R, U, Vs, l, gamma, has_tail = has_tail)

    return(ret)
}

#' Fit a GP with a mixed fidelity kernel.
#' @param lower_val How small may nonnegative quantities be?
#' @param smart_init If TRUE, initialize subspace using Principle Hessian Directions, lengthscales at mean squared distance, covars at 1 for main term, 0 for tail. Otherwise, randomly.
#' @param preinit Only useful is has_tail is TRUE. If pretail is TRUE, first fits a model without the tail and initializes on the result. If FALSE, standard init.
#' @return A list with optimization results. Note that optimization time includes preinitialization, if applicable.
gp_mf <- function(y, X, R, has_tail = TRUE, lower_val = 1e-4, smart_init = TRUE, preinit = FALSE) {
    N <- nrow(X)
    P <- ncol(X)

    # Some error catching
    if (R == 0 && !has_tail) {
        stop("Asking for no tail and a 0D active subspace yields an empty kernel.")
    }

    # Initialize 
    if (smart_init) {
        l_init <- rep(mean(dist(X)^2), P)
        if (R > 0) {
            u_init <- as.numeric(dr(y ~ X, method = 'phdy')$raw.evectors[,1:R])
        } else {
            u_init <- numeric(0)
        }
        init <- c(u_init, l_init)
        if (has_tail) {
            #gamma_init <- abs(rnorm(P-R+1,0,1e0))
            gamma_init <- c(var(y), rep(0, P-R+1))
        }  else {
            gamma_init <- var(y)
        }
        init <-  c(init, gamma_init)
    } else {
        l_init <- abs(rnorm(P))
        u_init <- rnorm(P*R)
        init <- c(u_init, l_init)
        if (has_tail) {
            gamma_init <- abs(rnorm(P-R+1,0,1e0))
        }  else {
            gamma_init <- abs(rnorm(1))
        }
        init <-  c(init, gamma_init)
    }

    if (preinit) {
        if (!has_tail) {
            warning("preinit only applicable to models with tails.")
        } else {
        ttp <- system.time(preoptim_ret <- optim(init, nll_wrapper, method = 'L-BFGS-B', 
              lower = c(rep(-Inf, R*P), rep(lower_val, P), 0, rep(0, has_tail*(P-R))), 
              upper = c(rep(Inf, R*P), rep(P, P), Inf, rep(Inf, has_tail*(P-R))), 
              y = y, X = X, R = R, has_tail = FALSE, control = list(trace = 0)))
        init <- preoptim_ret$par
        }
    }

    tt <- system.time(optim_ret <- optim(init, nll_wrapper, method = 'L-BFGS-B', 
          lower = c(rep(-Inf, R*P), rep(lower_val, P), 0, rep(0, has_tail*(P-R))), 
          upper = c(rep(Inf, R*P), rep(P, P), Inf, rep(Inf, has_tail*(P-R))), 
          y = y, X = X, R = R, has_tail = has_tail, control = list(trace = 0)))
    # Use trace = 3 for a good amount of output.

    # Create the return object.
    ret <- list(tt = tt, optim_ret = optim_ret)
    if (preinit) ret$tt <- ret$tt + ttp
    ret <- c(ret, par_2_param(optim_ret$par, P, R))
    ret$has_tail = has_tail
    ret$X <- X
    ret$y <- y
    ret$K <- form_cov(X, k = function(d) k(d, ret$U, ret$Vs, ret$l, ret$gamma, ret$has_tail), g = 1e-6)
    class(ret) <- "act_gp"

    return(ret)
}

#' 
predict.act_gp <- function(agp, XX) {
    K_xx_x <- form_cov(XX, agp$X, k = function(d) k(d, agp$U, agp$Vs, agp$l, agp$gamma, agp$has_tail), g = 1e-6)
    return(K_xx_x %*% solve(agp$K, agp$y))
}

print.act_gp <- function(agp) {
    cat(rep('-', as.integer(Sys.getenv("COLUMNS"))-1), '\n', sep = '')
    if (agp$has_tail) {
        cat("--Multi-Fidelity GP Kernel--\n")
    } else {
        cat("--Low-Rank GP Kernel--\n")
    }
    cat("Covariance Parameters:", agp$gamma, '\n')
    cat("Lengthscale Parameters:", agp$l, '\n')
    cat(rep('-', as.integer(Sys.getenv("COLUMNS"))-1), '\n', sep = '')
}
