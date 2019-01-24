#!/usr/bin/Rscript
#  R/kern_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.24.2019

#' The kernel to be studied
#' @param d The difference between two points.
#' @param U The high fidelity subspace, or, to be precise, a matrix the columns of which contain an orthonormal basis therefor.
#' @param Vs A list containing unit vectors filling the space after U.
#' @param l A numeric vector of nonnegative lengthscales. Gives the variance if we were to consider this a Gaussian.
#' @param gamma The covariance or scaling factors for each component of the sum; should be of length P - R + 1
k <- function(d, U, Vs, l, gamma) {
    #TODO: There are much nicer ways to do this.
    #Pu <- tcrossprod(U)
    R <- ncol(U)
    Pvs <- lapply(Vs, tcrossprod)
    corr <- gamma[1] * exp(-0.5 * t(d) %*% U %*% diag(1/l[1:R]) %*% t(U) %*% d)
    if (R < P) {
        for (p in (R+1):P) {
            corr <- corr + gamma[p-R+1] * exp(-0.5 * t(d) %*% Pvs[[p-R]] %*% d / l[p])
        }
    }
    return(corr)
}

#' Negative GP Log Likelihood.
#' @param y The response, a numeric.vector.
#' @param X The design, a matrix with as many rows as observations.
#' @param U The high fidelity subspace, or, to be precise, a matrix the columns of which contain an orthonormal basis therefor.
#' @param Vs A list containing unit vectors filling the space after U.
#' @param l A numeric vector of nonnegative lengthscales. Gives the variance if we were to consider this a Gaussian.
#' @param gamma The covariance or scaling factors for each component of the sum; should be of length P - R + 1
#' @param g The error variance.
nll <- function(y, X, U, Vs, l, gamma, g = 1e-6) {
    N <- nrow(X)
    P <- ncol(X)

    # Build the covariance matrix
    K <- matrix(NA, nrow = N, ncol = N)
    for (n1 in 1:N) {
        for (n2 in n1:N) {
            K[n1,n2] <- K[n2, n1] <- k(X[n1,] - X[n2,], U, Vs, l, gamma)
        }
    }
    diag(K) <- diag(K) + g

    # Evaluate loglik
    return(-dmvnorm(y, rep(0, N),  K, log = T))
}

#' A wrapper of the negative loglikelihood for use by optim
#' @param par is a vector, the first R*P elements give a matrix the range of which is a basis for the HF space. The next P parameters give lengthscales, and the P-R+1 parameters after that give covariances. As such, this vector should be of length (R*P+P+P-R+1). The first R*P are unconstrained, while the remaining should be bounded below by 0, or, preferably, a small constant.
nll_wrapper <- function(par, y, X, R) {
    N <- nrow(X)
    P <- ncol(X)

    # Form parameters
    L <- matrix(par[1:(R*P)], nrow = P)
    l <- par[(R*P+1):(R*P + P)]
    gamma <- par[(R*P+P+1):(R*P+P+P-R+1)]

    # Build orthonormal basis.
    B <- qr.Q(qr(cbind(L, matrix(rnorm((P-R)*P), nrow = P))))
    U <- B[,1:R]
    if (R < P) {
        Vs <- lapply((R+1):P, function(p) B[,p])
    } else {
        Vs <- NULL
    }

    return(nll(y, X, U, Vs, l, gamma))
}
