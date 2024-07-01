#### dependencies ####
library(gammaFrailty)
library(RcppClock)
library(dplyr)
library(tidyr)
library(purrr)
library(pbapply)
pbo = pboptions(type="txt")
# General folder
if(!dir.exists('gf/')){dir.create('gf/'); cat('[Create folder for Gamma Frailty results]') }

# Folder to store .Rda files
if(!dir.exists('gf/data_files')) {dir.create('gf/data_files'); cat('[Create data files folder]') }

# Folder to store plots
if(!dir.exists('gf/plots')) {dir.create('gf/plots'); cat('[Create plots folder]') }

# Run  with settingLab <- 'I' and p <- 20 for the smaller setting 
# Run  with settingLab <- 'II' and p <- 30 for the larger setting 

settingLab <- 'I'
if(!dir.exists(paste0('gf/data_files/', settingLab))){dir.create(paste0('gf/data_files/', settingLab)); cat('[Create setting-related data folder]')}

if(!dir.exists(paste0('gf/plots/', settingLab))){dir.create(paste0('gf/plots/', settingLab)); cat('[Create setting-related plots folder]')}


### Setting setup ####


# Setting parameters
p <- 20
q <- 8

xi <- 2/q
rho <- .5

m <- 0 # no covariates
set.seed(123)
int <- rep(c(-.25, .25), p/2)
b <-  rep(c(-.25,.25), m/2)
theta <- c(xi, rho, b, int)
repar_theta <- partorepar(theta)
d <- length(theta)
correlation_structure <- 'COMPOUND'
# save setting setup
save(p, q, xi, rho, m, int, b, theta, repar_theta, d, settingLab, correlation_structure, file = paste0('gf/data_files/', settingLab,'/setting_init.Rda'))
cat('[Setting setup saved]')


### Simulation setup ####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
set.seed(123)
par_init <- rep(0, length(theta))
sample_sizes <- c(2.5e3, 5e3, 10e3)
nsim <- 500
stepsizes <- c(0.5, 1, 2, 4, 8)

npasses <- 3
burn <- .25
range_lag <- p

strLab <- switch (correlation_structure,
                                                          'AR'     = 0,
                                                          'COMPOUND' = 1)
scls <- list(
     list('lab' = 'unitary', 'vec' = scls <- rep(1, d))
)


sim_settings <- expand_grid(
    sample_size = sample_sizes,
    stoc_seed = 1:nsim,
    stepsize = stepsizes,
    scaling_list = scls) %>%
    mutate(
        maxiter = sample_size * npasses,
        burn = round(sample_size * burn, 0),
        scaling = map_chr(scaling_list, ~.x$lab),
        scls = map(scaling_list, ~.x$vec)
    ) %>%
    select(-scaling_list)
save(sample_sizes, range_lag, par_init, nsim, stepsizes, sim_settings, npasses, burn, file = paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))
cat('[Simulation setup saved]')


ncores <- 1



#### ucminf ####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'ucminf'
path <- paste0('gf/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
    purrr::transpose(.l = sim_settings %>% select(-c(stepsize, scls, scaling, maxiter, burn)) %>% distinct()),
    function(x)try({

        n <- x$sample_size
        set.seed(n)
        X <- matrix(rbinom(m*n, 1, .5), n, m)

        data <- generate_data(
            INTERCEPT = int,
            BETA = b,
            X = X,
            Q = q,
            RHO = rho,
            SEED = x$stoc_seed,
            STRUCT = correlation_structure
        )

        suppressMessages(
            mod_obj <- fit_gammaFrailty2(
                DATA_LIST = list('DATA' = data, 'X' = X),
                METHOD = 'ucminf',
                CPP_CONTROL = list(),
                VERBOSEFLAG= 0,
                INIT = par_init,
                ITERATIONS_SUBSET = NULL,
                STRUCT = 'COMPOUND'
            )
        )

        mod_obj$lab <- lab
        mod_obj
}), cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))



### standard ####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'standard'
path <- paste0('gf/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
    purrr::transpose(.l = sim_settings ),
    function(x){

        n <- x$sample_size
        set.seed(n)
        X <- matrix(rbinom(m*n, 1, .5), n, m)

        data <- generate_data(
            INTERCEPT = int,
            BETA = b,
            X = X,
            Q = q,
            RHO = rho,
            SEED = x$stoc_seed,
            STRUCT = 'COMPOUND'
        )

        subtraj <- seq(0, npasses,.25)*n
        ctrl_sgd <- list(
            MAXT = x$maxiter,
            BURN = x$burn,
            STEPSIZE = x$stepsize,
            SCALEVEC = x$scls,
            PAR2 = x$stepsize,
            PAR3 = .501,
            NU = 1,
            STEPSIZEFLAG = 0,
            SAMPLING_WINDOW = 500,
            EACHCLOCK = 500
        )

        set.seed(x$stoc_seed)
        suppressMessages(
            mod_obj <- fit_gammaFrailty2(
                DATA_LIST = list('DATA' = data, 'X' = X),
                METHOD = 'standard',
                CPP_CONTROL = ctrl_sgd,
                PAIRS_RANGE = range_lag,
                VERBOSEFLAG= 0,
                INIT = par_init,
                ITERATIONS_SUBSET = subtraj,
                STRUCT = 'COMPOUND'
            )
        )


        mod_obj$lab <- lab
        mod_obj


        mod_obj$lab <- lab
        mod_obj
}, cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))

#### recycle_standard_500 ####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'recycle_standard_500'
path <- paste0('gf/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
    purrr::transpose(.l = sim_settings ),
    function(x){

        n <- x$sample_size
        set.seed(n)
        X <- matrix(rbinom(m*n, 1, .5), n, m)

        data <- generate_data(
            INTERCEPT = int,
            BETA = b,
            X = X,
            Q = q,
            RHO = rho,
            SEED = x$stoc_seed,
            STRUCT = 'COMPOUND'
        )

        subtraj <- seq(0, npasses,.25)*n
        ctrl_sgd <- list(
            MAXT = x$maxiter,
            BURN = x$burn,
            STEPSIZE = x$stepsize,
            SCALEVEC = x$scls,
            PAR2 = x$stepsize,
            PAR3 = .501,
            NU = 1,
            STEPSIZEFLAG = 0,
            SAMPLING_WINDOW = 500,
            EACHCLOCK = 500,
            SEED = x$stoc_seed
        )

        set.seed(x$stoc_seed)
        suppressMessages(
            mod_obj <- fit_gammaFrailty2(
                DATA_LIST = list('DATA' = data, 'X' = X),
                METHOD = 'recycle_standard',
                CPP_CONTROL = ctrl_sgd,
                PAIRS_RANGE = range_lag,
                VERBOSEFLAG= 0,
                INIT = par_init,
                ITERATIONS_SUBSET = subtraj,
                STRUCT = 'COMPOUND'
            )
        )


        mod_obj$lab <- lab
        mod_obj
    }, cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))






#### recycle_hyper_100 ####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'recycle_hyper'
path <- paste0('gf/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
    purrr::transpose(.l = sim_settings ),
    function(x){

        n <- x$sample_size
        set.seed(n)
        X <- matrix(rbinom(m*n, 1, .5), n, m)

        data <- generate_data(
            INTERCEPT = int,
            BETA = b,
            X = X,
            Q = q,
            RHO = rho,
            SEED = x$stoc_seed,
            STRUCT = 'COMPOUND'
        )

        subtraj <- seq(0, npasses,.25)*n
        ctrl_sgd <- list(
            MAXT = x$maxiter,
            BURN = x$burn,
            STEPSIZE = x$stepsize,
            SCALEVEC = x$scls,
            PAR2 = x$stepsize,
            PAR3 = .501,
            NU = 1,
            STEPSIZEFLAG = 0,
            SAMPLING_WINDOW = 100,
            EACHCLOCK = 500
        )

        set.seed(x$stoc_seed)
        suppressMessages(
            mod_obj <- fit_gammaFrailty2(
                DATA_LIST = list('DATA' = data, 'X' = X),
                METHOD = 'recycle_hyper',
                CPP_CONTROL = ctrl_sgd,
                PAIRS_RANGE = range_lag,
                VERBOSEFLAG= 0,
                INIT = par_init,
                ITERATIONS_SUBSET = subtraj,
                STRUCT = 'COMPOUND'
            )
        )


        mod_obj$lab <- lab
        mod_obj
    }, cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))





#### recycle_hyper_500 ####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'recycle_hyper_500'
path <- paste0('gf/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
    purrr::transpose(.l = sim_settings ),
    function(x){

        n <- x$sample_size
        set.seed(n)
        X <- matrix(rbinom(m*n, 1, .5), n, m)

        data <- generate_data(
            INTERCEPT = int,
            BETA = b,
            X = X,
            Q = q,
            RHO = rho,
            SEED = x$stoc_seed,
            STRUCT = 'COMPOUND'
        )

        subtraj <- seq(0, npasses,.25)*n
        ctrl_sgd <- list(
            MAXT = x$maxiter,
            BURN = x$burn,
            STEPSIZE = x$stepsize,
            SCALEVEC = x$scls,
            PAR2 = x$stepsize,
            PAR3 = .501,
            NU = 1,
            STEPSIZEFLAG = 0,
            SAMPLING_WINDOW = 500,
            EACHCLOCK = 500,
            SEED = x$stoc_seed
        )

        set.seed(x$stoc_seed)
        suppressMessages(
            mod_obj <- fit_gammaFrailty2(
                DATA_LIST = list('DATA' = data, 'X' = X),
                METHOD = 'recycle_hyper',
                CPP_CONTROL = ctrl_sgd,
                PAIRS_RANGE = range_lag,
                VERBOSEFLAG= 0,
                INIT = par_init,
                ITERATIONS_SUBSET = subtraj,
                STRUCT = 'COMPOUND'
            )
        )


        mod_obj$lab <- lab
        mod_obj
    }, cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))

#### recycle_hyper_1000 ####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'recycle_hyper_1000'
path <- paste0('gf/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
    purrr::transpose(.l = sim_settings ),
    function(x){

        n <- x$sample_size
        set.seed(n)
        X <- matrix(rbinom(m*n, 1, .5), n, m)

        data <- generate_data(
            INTERCEPT = int,
            BETA = b,
            X = X,
            Q = q,
            RHO = rho,
            SEED = x$stoc_seed,
            STRUCT = 'COMPOUND'
        )

        subtraj <- seq(0, npasses,.25)*n
        ctrl_sgd <- list(
            MAXT = x$maxiter,
            BURN = x$burn,
            STEPSIZE = x$stepsize,
            SCALEVEC = x$scls,
            PAR2 = x$stepsize,
            PAR3 = .501,
            NU = 1,
            STEPSIZEFLAG = 0,
            SAMPLING_WINDOW = 1000,
            EACHCLOCK = 500,
            SEED = x$stoc_seed
        )

        set.seed(x$stoc_seed)
        suppressMessages(
            mod_obj <- fit_gammaFrailty2(
                DATA_LIST = list('DATA' = data, 'X' = X),
                METHOD = 'recycle_hyper',
                CPP_CONTROL = ctrl_sgd,
                PAIRS_RANGE = range_lag,
                VERBOSEFLAG= 0,
                INIT = par_init,
                ITERATIONS_SUBSET = subtraj,
                STRUCT = 'COMPOUND'
            )
        )


        mod_obj$lab <- lab
        mod_obj
    }, cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))




#### COLLECT #####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))

true_tib <- tibble(par = 1:length(repar_theta), true_val = repar_theta) #%>%  mutate_at(vars(par), as.factor)
load(paste0('gf/data_files/', settingLab, '/ucminf/est.Rda'))
num_est <- sim_settings %>%
    mutate(stepsize = NA, scaling = NA, scls = NA, maxiter = NA, burn = NA, meth = 'ucminf') %>%
    distinct() %>%
    mutate(mod = est_obj_list)

load(paste0('gf/data_files/', settingLab, '/recycle_standard_500/est.Rda'))
recycle_standard_500_est <- sim_settings %>%
    mutate(meth = 'recycle_standard_500') %>%
    distinct() %>%
    mutate(mod = est_obj_list)

load(paste0('gf/data_files/', settingLab, '/recycle_hyper/est.Rda'))
recycle_hyper_est <- sim_settings %>%
    mutate(meth = 'recycle_hyper_100') %>%
    distinct() %>%
    mutate(mod = est_obj_list)

load(paste0('gf/data_files/', settingLab, '/recycle_hyper_500/est.Rda'))
recycle_hyper_500_est <- sim_settings %>%
    mutate(meth = 'recycle_hyper_500') %>%
    distinct() %>%
    mutate(mod = est_obj_list)
load(paste0('gf/data_files/', settingLab, '/recycle_hyper_1000/est.Rda'))
recycle_hyper_1000_est <- sim_settings %>%
    mutate(meth = 'recycle_hyper_1000') %>%
    distinct() %>%
    mutate(mod = est_obj_list)


stoc_est <- recycle_standard_500_est %>%
    bind_rows(recycle_hyper_est) %>%
    bind_rows(recycle_hyper_500_est) %>%
    bind_rows(recycle_hyper_1000_est)

path_tab <- stoc_est %>%
    mutate( err= map_lgl(mod, ~is.null(.x))) %>%
    filter(err == F) %>%
    mutate(path_theta = map(mod, ~get_tidy_path(.x, 'path_av_theta', F))) %>%
    select(-mod)

save(num_est, stoc_est, path_tab, true_tib, file = paste0('gf/data_files/', settingLab, '/est.Rda'))

#### variance ####
load(paste0('gf/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf/data_files/',settingLab,'/sim_setup.Rda'))
load(paste0('gf/data_files/',settingLab,'/est.Rda'))
selected_passes <- seq(.5,npasses,.5)

#### Computation ####
path_selected <- path_tab %>%
    unnest(path_theta) %>%
    mutate(pass = iter/sample_size) %>%
    filter(pass %in% selected_passes) %>%
    distinct()

var_list <- pbapply::pblapply(path_selected %>% purrr::transpose(),
    function(x){

        # regenerate data
        n <- x$sample_size
        set.seed(n)
        X <- matrix(rbinom(m*n, 1, .5), n, m)
        data <- generate_data(
            INTERCEPT = int,
            BETA = b,
            X = X,
            Q = q,
            RHO = rho,
            SEED = x$stoc_seed,
            STRUCT = correlation_structure
        )

        # read method
        method <- switch (x$meth,
                          'GD'     = 0,
                          'ucminf' = 0,
                          'recycle_standard_500'= 1,
                          'recycle_hyper_100'   = 2,
                          'recycle_hyper_500'   = 2,
                          'recycle_hyper_1000'  = 2)
        strLab <- switch (correlation_structure,
                          'AR'     = 0,
                          'COMPOUND' = 1)

        out <- sampleVar(
            THETA = x$path_av_theta,
            DATA = data,
            X = X,
            NU = 1,
            METHOD = method,
            RANGE = x$iter - x$burn,
            TOTFLAG = T,
            PRINTFLAG = F,
            STRUCT = strLab

        )

        se <- out$se


        vt_tot <- out$var$var_tot %>% diag() %>% sum()
        vt_stoc <- out$var$var_stoc %>% diag() %>% sum()
        vt_stat <- out$var$var_stoc %>% diag() %>% sum()
        tr <- tibble( var_type = c('stoc', 'stat', 'tot'), tr = c(vt_stoc, vt_stat, vt_tot))


        return(list('se'=se, 'trace' = tr ))
    },
    cl = ncores)

path_selected$asy_var <- var_list
save(path_selected, file = paste0('gf/data_files/', settingLab,'/','path_variance.Rda'))


var_list_num <- pbapply::pblapply(num_est %>% purrr::transpose(),
                                  function(x){

                                      # regenerate data
                                      n <- x$sample_size
                                      set.seed(n)
                                      X <- matrix(rbinom(m*n, 1, .5), n, m)
                                      data <- generate_data(
                                          INTERCEPT = int,
                                          BETA = b,
                                          X = X,
                                          Q = q,
                                          RHO = rho,
                                          SEED = x$stoc_seed,
            STRUCT = correlation_structure
                                      )

                                      # read method
                                      method <- switch (x$meth,
                                                        'GD'     = 0,
                                                        'ucminf' = 0,
                                                        'recycle_standard_500'= 1,
                                                        'recycle_hyper_100'   = 2,
                                                        'recycle_hyper_500'   = 2,
                                                        'recycle_hyper_1000'  = 2)
                                      strLab <- switch (correlation_structure,
                                                          'AR'     = 0,
                                                          'COMPOUND' = 1)

                                      out <- sampleVar(
                                          THETA = x$mod$fit$par,
                                          DATA = data,
                                          X = X,
                                          NU = 1,
                                          METHOD = method,
                                          RANGE = 1,
                                          TOTFLAG = T,
                                          PRINTFLAG = F,
                                          STRUCT = strLab
                                      )

                                      se <- out$se


                                      vt_tot <- out$var$var_tot %>% diag() %>% sum()
                                      vt_stoc <- out$var$var_stoc %>% diag() %>% sum()
                                      vt_stat <- out$var$var_stoc %>% diag() %>% sum()
                                      tr <- tibble( var_type = c('stoc', 'stat', 'tot'), tr = c(vt_stoc, vt_stat, vt_tot))


                                      return(list('se'=se, 'trace' = tr ))
                                  },
                                  cl = ncores)
num_est$asy_var <- var_list_num

save(num_est, file = paste0('gf/data_files/', settingLab,'/','num_variance.Rda'))







