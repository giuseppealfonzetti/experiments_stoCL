% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sampleH}
\alias{sampleH}
\title{Compute sample covariance matrix}
\usage{
sampleH(THETA, DATA, CONSTRAINTS, INVERTFLAG = FALSE, VERBOSEFLAG = FALSE)
}
\arguments{
\item{THETA}{Parameter vector}

\item{DATA}{Matrix with \code{n} rows and \code{p} columns}

\item{CONSTRAINTS}{Vector of booleans: TRUE denotes free-to-estimate parameters}

\item{INVERTFLAG}{TRUE for inverse of the sample Hessian}

\item{VERBOSEFLAG}{Verbose output}
}
\description{
Compute sample covariance matrix
}
