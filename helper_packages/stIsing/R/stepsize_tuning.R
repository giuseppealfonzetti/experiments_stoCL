
#'@export
stepsize_tuning4 <- function(
        DATA_LIST = list('DATA', 'CONSTRAINTS', 'HOLDOUT'),
        METHOD,
        CPP_CONTROL = list(),
        STEPSIZE_INIT = NULL,
        LENGTH = 0.5,
        INIT = NULL,
        VERBOSEFLAG = 0,
        MAXATTEMPT = 10
){
    if(require(dplyr)){
        out <- list()
        start_time <- Sys.time()

        p <- ncol(DATA_LIST$DATA)
        d <- p + p*(p-1)/2

        # Check Initialisation
        if(is.vector(INIT)){
            if(length(INIT)!=d)
                stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
            else
                #message('1. Initialising at init vector.')
            out$theta_init <-  INIT
        }else{
            if(is.null(INIT))
                #message('1. Initialising at zero vector.')
            out$theta_init <-  rep(0, d)
        }

        # Check if method entry is correct
        if(!(METHOD %in% c('standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper'))) stop('"METHOD" not available.')
        out$method <- METHOD

        # Check holdout
        if(is.null(DATA_LIST$HOLDOUT)){
            CPP_CONTROL$HOLDOUTFLAG <- F
            DATA_LIST$HOLDOUT <- matrix(0, 1, p)
        }else{
            if(ncol(DATA_LIST$HOLDOUT)!=ncol(DATA_LIST$DATA)) stop(' "DATA" and "HOLDOUT" must have the same number of nodes.')
            CPP_CONTROL$HOLDOUTFLAG <- T

        }

        # Check stochastic control parameters
        cpp_ctrl <- check_stoc_args(CPP_CONTROL, N = n, D = d)

        # Guarantee reproducibility stochastic optimisation with bernoulli sampling
        # For the other schemes the seed is passed directly to the cpp function
        set.seed(cpp_ctrl$SEED)

        # Collect and rearrange arguments to pass to cpp function
        args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

        if(METHOD == 'standard'){args$METHODFLAG <- 1} else if(METHOD == 'bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'hyper'){args$METHODFLAG <- 3}else if(METHOD == 'recycle_standard'){args$METHODFLAG <- 4}else if(METHOD == 'recycle_hyper'){args$METHODFLAG <- 5}

        # objective function on holdout data
        Rwr_ncl <- function(par){ ncl(DATA_LIST$HOLDOUT, par, DATA_LIST$CONSTRAINTS)$nll/nrow(DATA_LIST$HOLDOUT) }

        cond <- T
        args$STEPSIZE <- STEPSIZE_INIT
        fit <- do.call(isingGraph3, args)
        nll <- Rwr_ncl(fit$path_av_theta[[length(fit$path_av_theta)]])
        tib <- tibble(stepsize = args$STEPSIZE, hnll = nll)
        tib$theta <- list(fit$path_av_theta[[length(fit$path_av_theta)]])
        tib$path_theta <- list(tibble(iter = fit$iter_idx, path_av_theta = fit$path_av_theta))

        res <- tib
        new_args <- args
        for (i in 1:MAXATTEMPT) {
            new_args$STEPSIZE <- new_args$STEPSIZE*LENGTH

            new_fit <- do.call(isingGraph3, new_args)
            new_nll <- Rwr_ncl(new_fit$path_av_theta[[length(new_fit$path_av_theta)]])


            tib <- tibble(stepsize = new_args$STEPSIZE,  hnll = new_nll)
            tib$theta <- list(new_fit$path_av_theta[[length(new_fit$path_av_theta)]])
            tib$path_theta <- list(tibble(iter = new_fit$iter_idx, path_av_theta = new_fit$path_av_theta))
            res <- res %>% bind_rows(tib)

            cond <- (new_nll < nll)

            if(!is.logical(cond)){stop('likelihood differene undefined')}
            if(!cond) break;
            nll <- new_nll

        }



        best = res |> filter(hnll == min(hnll))

        return(list(
            best_theta = best$theta[[1]],
            grid = res,
            best_step = best$stepsize[[1]],
            best_nll = best$hnll[[1]],
            best_path = best$path_theta[[1]]
        )
        )
    }else{break}

}
