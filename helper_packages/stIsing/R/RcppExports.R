# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Negative composite log-likelihood
#'
#' Compute the negative composite log-likelihood of an Ising graph and its gradient
#'
#' @param DATA Data matrix with `n` rows and `p` columns
#' @param THETA Parameter vector.
#' @param CONSTRAINTS Vector of booleans: TRUE denotes free-to-estimate parameters
#' @param VERBOSEFLAG Verbose output.
#'
#' @export
ncl <- function(DATA, THETA, CONSTRAINTS, VERBOSEFLAG = FALSE) {
    .Call(`_stIsing_ncl`, DATA, THETA, CONSTRAINTS, VERBOSEFLAG)
}

isingGraph3 <- function(DATA, HOLDOUT, THETA_INIT, CONSTRAINTS, MAXT, BURN, STEPSIZE, NU, METHODFLAG, SCALEVEC, SEED = 123L, VERBOSEFLAG = FALSE, HOLDOUTFLAG = FALSE, PAR1 = 1, PAR2 = 1, PAR3 = .501, SAMPLING_WINDOW = 1L, EACH = 1L, EACHCLOCK = 100L, STEPSIZEFLAG = 0L, T_INIT = 1L) {
    .Call(`_stIsing_isingGraph3`, DATA, HOLDOUT, THETA_INIT, CONSTRAINTS, MAXT, BURN, STEPSIZE, NU, METHODFLAG, SCALEVEC, SEED, VERBOSEFLAG, HOLDOUTFLAG, PAR1, PAR2, PAR3, SAMPLING_WINDOW, EACH, EACHCLOCK, STEPSIZEFLAG, T_INIT)
}

rmultinom_wrapper <- function(prob, classes, batch, K) {
    .Call(`_stIsing_rmultinom_wrapper`, prob, classes, batch, K)
}

hyper_sampling <- function(K, N, SEED) {
    .Call(`_stIsing_hyper_sampling`, K, N, SEED)
}

unit_sampling <- function(N, SEED) {
    .Call(`_stIsing_unit_sampling`, N, SEED)
}

components_given_unit <- function(UNIT, K) {
    .Call(`_stIsing_components_given_unit`, UNIT, K)
}

bernoulli_sampling <- function(K, N, PROB) {
    .Call(`_stIsing_bernoulli_sampling`, K, N, PROB)
}

index_to_component <- function(P, N, INDEX) {
    .Call(`_stIsing_index_to_component`, P, N, INDEX)
}

#' Compute sample covariance matrix
#'
#' @param THETA Parameter vector
#' @param DATA Matrix with `n` rows and `p` columns
#' @param CONSTRAINTS Vector of booleans: TRUE denotes free-to-estimate parameters
#' @param INVERTFLAG TRUE for inverse of the sample Hessian
#' @param VERBOSEFLAG Verbose output
#'
#' @export
sampleH <- function(THETA, DATA, CONSTRAINTS, INVERTFLAG = FALSE, VERBOSEFLAG = FALSE) {
    .Call(`_stIsing_sampleH`, THETA, DATA, CONSTRAINTS, INVERTFLAG, VERBOSEFLAG)
}

#' Compute sample variability matrix
#'
#' @param THETA Parameter vector
#' @param DATA Matrix with `n` rows and `p` columns
#' @param CONSTRAINTS Vector of booleans: TRUE denotes free-to-estimate parameters
#' @param VERBOSEFLAG Verbose output
#'
#' @export
sampleJ <- function(THETA, DATA, CONSTRAINTS, VERBOSEFLAG = FALSE) {
    .Call(`_stIsing_sampleJ`, THETA, DATA, CONSTRAINTS, VERBOSEFLAG)
}

#' Compute sample covariance matrix
#'
#' @param THETA Parameter vector
#' @param DATA Matrix with `n` rows and `p` columns
#' @param CONSTRAINTS Vector of booleans: TRUE denotes free-to-estimate parameters
#' @param NU Number of pairs per iteration on average.
#' @param METHOD Allowed choices are "ucminf", "standard", "bernoulli", "hyper",
#' "recycle_standard", "recycle_bernoulli", "recycle_hyper"
#' @param RANGE Number of iterations after the burn-in.
#' @param TOTFLAG Compute total variance
#' @param PRINTFLAG Verbose output
#'
#' @export
sampleVar <- function(THETA, DATA, CONSTRAINTS, NU, METHOD, RANGE, TOTFLAG, PRINTFLAG) {
    .Call(`_stIsing_sampleVar`, THETA, DATA, CONSTRAINTS, NU, METHOD, RANGE, TOTFLAG, PRINTFLAG)
}

