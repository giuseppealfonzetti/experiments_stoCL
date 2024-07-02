generate_C <- function(RHO, P, STRUCT){
    if(STRUCT=='AR'){
        out <- diag(rep(1,P))
        for (j in 2:P) {
            for (k in 1:(j-1)) {
                out[j,k] <- RHO^(abs(j-k)/2)
            }
        }
        out <- out+t(out)-diag(rep(1,P))
    }else if(STRUCT=='COMPOUND'){
        out <- matrix(rep(RHO^.5, P, P), P, P);
        diag(out) <- 1
    }
    return(out)
}

generate_mgamma <- function(Q, C, SEED){
    set.seed(SEED)
    out <- colMeans((mvtnorm::rmvnorm(n = Q, sigma = C))^2)

    return(out)
}

#' Simulate a synthetic dataset
#'
#' @param INTERCEPT Intercept vector of length `p`
#' @param BETA Covariates coefficients vector
#' @param X External covariates matrix with `n` rows
#' @param Q Value of q in the data generating process. See Henderson and Shimakura (2003)
#' @param RHO Correlation
#' @param STRUCT Structure of the correlation matrix. Allowed values are "AR" or "COMPOUND".
#' @param SEED Seed for reproducibility.
#'
#'@export
generate_data <- function(INTERCEPT, BETA, X, Q, RHO, SEED, STRUCT){

    if(!is.matrix(X)){cat('Error: X must be a matrix!\n'); return(NULL)}
    if(length(BETA)!=ncol(X)){cat('Error: Dimensions of beta and X not compatible!\n'); return(NULL)}
    if(RHO < 0){cat('Error: rho must be non-negative!\n'); return(NULL)}

    p <- length(INTERCEPT)
    C <- generate_C(RHO = RHO, P = p, STRUCT = STRUCT)
    n <- nrow(X)
    Z <- purrr::reduce(purrr::map(1:n, ~generate_mgamma(Q, C, SEED = SEED + .x)), rbind); rownames(Z) <- NULL
    u <- t(sapply(1:n, function(i) exp(INTERCEPT+ as.numeric(crossprod(BETA, X[i,])))))

    set.seed(SEED)
    Zu <- Z*u

    out <- apply(Zu, c(1, 2), function(x)rpois(1, lambda = x))

    return(out)
}

#' #'@export
#' uni_dnbinom <- function(u, n, eps){
#'     nf <- factorial(n)
#'
#'     pr <- 1
#'     for (i in 0:(n-1)) {
#'         pr <- pr * (1+i*eps)
#'     }
#'
#'     frac <- (u^n)/(nf *(1+eps*u)^(n+1/eps))
#'     out <- frac * pr
#'
#'     return(as.numeric(out))
#' }
