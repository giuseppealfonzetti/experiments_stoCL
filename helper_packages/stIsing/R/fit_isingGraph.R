utils::globalVariables(c("clock"))
#' Ising model estimation
#'
#' @param DATA_LIST List containing:
#' \itemize{
#'  \item{`DATA`}{ Data matrix with `n` rows and `p` columns}
#'  \item{`CONSTRAINTS`}{ Vector of booleans: TRUE denotes free-to-estimate parameters }
#' }
#' @param METHOD Allowed choices are "ucminf", "standard", "bernoulli", "hyper",
#' "recycle_standard", "recycle_bernoulli", "recycle_hyper".
#' @param CPP_CONTROL List of arguments to be passed to the stochastic optimizer:
#' \itemize{
#'  \item{`MAXT`}{ Number of iterations.}
#'  \item{`BURN`}{ Scalar between 0 and 1 denoting the share of `n` to be used as burn-in iterations.}
#'  \item{`STEPSIZE`}{ Initial stepsize parameter.}
#'  \item{`NU`}{ Number of pairs per iteration on average.}
#'  \item{`PAR1`}{ Hyperparameter for stepsize scheduling by Xu (2011): Scaling.}
#'  \item{`PAR2`}{ Hyperparameter for stepsize scheduling by Xu (2011): Smallest Hessian eigenvalue.}
#'  \item{`PAR3`}{ Hyperparameter for stepsize scheduling by Xu (2011): Decay rate}
#'  \item{`STEPSIZEFLAG`}{ Choose stepsize scheduling: Set 0 for Polyak and Juditsky (1992), 1 for Xu (2011).}
#'  \item{`SAMPLING_WINDOW`}{ How many iterations to wait before resampling}
#'  \item{`EACHCLOCK`}{ How often (in terms of iterations) to measure single iteration computational times (using \code{RcppClock})}
#'  \item{`SEED`}{ Random seed for reproducibility.}
#' }
#' @param UCMINF_CONTROL List of arguments to be passed to \link[ucminf]{ucminf}:
#' \itemize{
#'  \item{"ctrl"}{ Passed to `ctrl` argument in \link[ucminf]{ucminf}}
#'  \item{"hessian"}{ Passed to `hessian` argument in \link[ucminf]{ucminf}}
#' }
#' @param INIT Initialization vector.
#' @param VERBOSEFLAG Verbose output.
#'
#' @export
fit_isingGraph3 <- function(
        DATA_LIST = list('DATA', 'CONSTRAINTS', 'HOLDOUT'),
        METHOD,
        CPP_CONTROL = list(
            MAXT = 1000,
            BURN = 500,
            STEPSIZE = .01,
            STEPSIZE0 = NULL,
            NU = 1,
            SEED = 123
        ),
        UCMINF_CONTROL = list('ctrl' = list(), 'hessian' = 0 ),
        INIT = NULL,
        VERBOSEFLAG = 0
){

    out <- list()
    start_time <- Sys.time()

    p <- ncol(DATA_LIST$DATA)
    d <- p + p*(p-1)/2
    n <- nrow(DATA_LIST$DATA)

    # Check Initialisation
    if(is.vector(INIT)){
        if(length(INIT)!=d)
            stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
        else
            message('1. Initialising at init vector.')
        out$theta_init <-  INIT
    }else{
        if(is.null(INIT))
            message('1. Initialising at zero vector.')
        out$theta_init <-  rep(0, d)
    }

    # Check if method entry is correct
    if(!(METHOD %in% c('ucminf','GD', 'standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper'))) stop('"METHOD" not available.')
    out$method <- METHOD

    # Check holdout
    if(is.null(DATA_LIST$HOLDOUT)){
        CPP_CONTROL$HOLDOUTFLAG <- F
        DATA_LIST$HOLDOUT <- matrix(0, 1, p)
    }else{
        if(ncol(DATA_LIST$HOLDOUT)!=ncol(DATA_LIST$DATA)) stop(' "DATA" and "HOLDOUT" must have the same number of nodes.')
        CPP_CONTROL$HOLDOUTFLAG <- T

    }

    # Numerical optimisation
    if(METHOD == 'ucminf'){

        message('2. Optimising with ucminf...')
        # R wrapper of cpp function for negative composite log-likelihood
        Rwr_ncl <- function(par){ ncl(DATA_LIST$DATA, par, DATA_LIST$CONSTRAINTS)$nll }
        #as.matrix(data), par, Q, VERBOSEFLAG = F
        # R wrapper of cpp function for negative composite score
        Rwr_ngr <- function(par){ ncl(DATA_LIST$DATA, par, DATA_LIST$CONSTRAINTS)$ngradient }

        # list of ucminf args
        args <- list(
            'par' = out$theta_init,
            'fn' = Rwr_ncl,
            'gr' = Rwr_ngr,
            'control' = UCMINF_CONTROL$ctrl,
            'hessian' = UCMINF_CONTROL$hessian)

        # optimisation
        opt <- do.call(ucminf::ucminf, args)
        out$fit <- opt

        out$control <- UCMINF_CONTROL
        out$theta   <- opt$par

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
        message('3. Done! (', round(out$time,2),' secs)')
        return(out)
    }

    # Stochastic approximation of numerical optimiser
    if(METHOD %in% c('standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper')){

        message(paste0('2. Optimising with ', METHOD, '...'))

        # Check stochastic control parameters
        cpp_ctrl <- check_stoc_args(CPP_CONTROL, N = n, D = d)

          set.seed(cpp_ctrl$SEED)

        # Collect and rearrange arguments to pass to cpp function
        args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )


        if(METHOD == 'standard'){args$METHODFLAG <- 1} else if(METHOD == 'bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'hyper'){args$METHODFLAG <- 3}else if(METHOD == 'recycle_standard'){args$METHODFLAG <- 4}else if(METHOD == 'recycle_hyper'){args$METHODFLAG <- 5}

        fit <- do.call(isingGraph3, args)

        message(paste0('3. Storing results'))

        out$theta <- fit$path_av_theta[[length(fit$path_av_theta)]]




        fit$methodflag <- NULL



        out$control <- cpp_ctrl
        out$fit <- fit
        if(requireNamespace('RcppClock', quietly = TRUE)) out$clock <- summary(clock, units = 's')

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

        message('4. Done! (', round(out$time,2),' secs)')
        return(out)

    }
}
