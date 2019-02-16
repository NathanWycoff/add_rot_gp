/** src/kern_funcs.cpp Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.26.2019  */

#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
double fast_k_cpp(NumericVector d, NumericVector gamma, int R, bool has_tail) {
    double corr = 0;
    int P = d.size();

    // Get sum of active subspace
    double as_sum = 0;
    for (int r = 0; r < R; r++) {
        as_sum += d(r)*d(r);
    }

    if (R > 0) {
        corr += gamma(0) * exp(-0.5 * as_sum);
    }
    if (has_tail && R < P) {
        for (int p = R; p < P; p++) {
            corr += gamma(p-R+1) * exp(-0.5 * d(p)*d(p));
        }
    }

    return(corr);
}

// [[Rcpp::export]]
NumericMatrix form_cov_cpp(NumericMatrix X, NumericVector gamma, int R, bool has_tail, double g) {
    NumericMatrix K(X.nrow(), X.nrow());
    for (int n1 = 0; n1 < X.nrow(); n1++) {
        for (int n2 = n1; n2 < X.nrow(); n2++) {
            K(n1,n2) = K(n2,n1) = fast_k_cpp(X(n1,_) - X(n2,_), gamma, R, has_tail);
            if (n1 == n2) {
                K(n1,n1) += g;
            }
        }
    }
    return(K);
}

// [[Rcpp::export]]
NumericMatrix form_crosscov_cpp(NumericMatrix X1, NumericMatrix X2, NumericVector gamma, int R, bool has_tail, double g) {
    NumericMatrix K(X1.nrow(), X2.nrow());
    for (int n1 = 0; n1 < X1.nrow(); n1++) {
        for (int n2 = 0; n2 < X2.nrow(); n2++) {
            K(n1,n2) = fast_k_cpp(X1(n1,_) - X2(n2,_), gamma, R, has_tail);
        }
    }
    return(K);
}


//#' Form GP Covariance structure. Nugget added to diagonal if X2 not given.
//form_cov <- function(X1, X2 = NULL, k, g = 1e-6) {
//    if (is.null(X2)) {
//        K <- matrix(NA, nrow = nrow(X1), ncol = nrow(X1))
//        for (n1 in 1:nrow(X1)) {
//            for (n2 in n1:nrow(X1)) {
//                K[n1,n2] <- K[n2,n1] <- k(X1[n1,] - X1[n2,])
//            }
//        }
//        diag(K) <- diag(K) + g
//    } else {
//        K <- matrix(NA, nrow = nrow(X1), ncol = nrow(X2))
//        for (n1 in 1:nrow(X1)) {
//            for (n2 in 1:nrow(X2)) {
//                K[n1,n2] <- k(X1[n1,] - X2[n2,])
//            }
//        }
//    }
//    return(K)
//}
