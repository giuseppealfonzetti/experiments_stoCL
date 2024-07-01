#### dependencies ###########
library(stIsing)
library(RcppClock)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)
library(IsingSampler)

# General folder
if(!dir.exists('isi/')){dir.create('isi/'); cat('[Create folder for Ising results]') }

# Folder to store .Rda files
if(!dir.exists('isi/data_files')) {dir.create('isi/data_files'); cat('[Create data files folder]') }

# Folder to store plots
if(!dir.exists('isi/plots')) {dir.create('isi/plots'); cat('[Create plots folder]') }


## Specify setting label:

# Run  with settingLab <- 'I' and p <- 10 for the smaller setting 
# Run  with settingLab <- 'II' and p <- 20 for the larger setting 

settingLab <- 'I'
if(!dir.exists(paste0('isi/data_files/', settingLab))){dir.create(paste0('isi/data_files/', settingLab)); cat('[Create setting-related data folder]')}

if(!dir.exists(paste0('isi/plots/', settingLab))){dir.create(paste0('isi/plots/', settingLab)); cat('[Create setting-related plots folder]')}


### Setting setup ####


# Setting parameters
p <- 10
d <- p + p*(p-1)/2

parNodes <- rep(c(-.5,.5), p/2)
parEdgesMat <- matrix(0,p, p)
counter <- 1
for (col in 1:(p-1)) {
    for (row in (col+1):p) {
        if(abs(col-row)==1 & (min(col,row)!=p/2)) parEdgesMat[row, col] <- .5
        if(abs(col-row)==p/2) parEdgesMat[row, col] <- -.5
    }
}

true_theta <- parNodes
for (col in 1:(p-1)) {
    for (row in (col+1):p) {
        true_theta <- c(true_theta, parEdgesMat[row, col])
    }
}
true_theta
Q <- rep(TRUE, length(true_theta))
graph_mat <- ising_from_theta_to_emat(true_theta, p)
graph_mat <- graph_mat + t(graph_mat)
thr_vec <- true_theta[1:p]

# save setting setup
save(p, d, parNodes, parEdgesMat, true_theta, Q, graph_mat, thr_vec, settingLab,  file = paste0('isi/data_files/', settingLab,'/setting_init.Rda'))
cat('[Setting setup saved]')

### Simulation setup ####

load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
set.seed(123)
par_init <- rep(0, d)

sample_sizes <- c(2.5e3, 5e3, 10e3)
nsim <- 500
stepsizes <- c(.25, .5, 1, 2, 4)
npasses <- 3
burn <- .25

set.seed(min(sample_sizes))
data <- IsingSampler::IsingSampler(min(sample_sizes), graph_mat, thr_vec, 1, method = "direct")

# Option for additional multidimensional stepsize scaling. Not used in the paper so set to unit vector
scls <- list(
    list('lab' = 'unitary', 'vec' = scls <- rep(1, d))
)

#'

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
save(sample_sizes, par_init, nsim, stepsizes, sim_settings, npasses, burn, file = paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))
cat('[Simulation setup saved]')


# select the number of cores
ncores <- 1



#### ucminf ####
load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'ucminf'
path <- paste0('isi/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
  purrr::transpose(.l = sim_settings %>% select(-c(stepsize, scls, scaling, maxiter, burn)) %>% distinct()),
  function(x)try({
    
    n <- x$sample_size
    set.seed(x$stoc_seed)
    data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
    
    
    suppressMessages(
      mod_obj <- fit_isingGraph3(
        DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
        METHOD = 'ucminf',
        CPP_CONTROL = list(),
        INIT = par_init,
        VERBOSEFLAG = 0
      )
    )
    
    mod_obj$lab <- lab
    mod_obj
  }), cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))



#### STANDARD ####
load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'standard'
path <- paste0('isi/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
  purrr::transpose(.l = sim_settings ),
  function(x)try({
    
    n <- x$sample_size
    set.seed(x$stoc_seed)
    data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
    
    subtraj <- seq(0, x$maxiter, 250)
    cpp_ctrl <- list(
      MAXT = x$maxiter,
      BURN = x$burn,
      STEPSIZE = x$stepsize,
      SCALEVEC = x$scls,
      PAR1 = 1,
      PAR2 = x$stepsize,
      PAR3 = .501,
      NU = 1,
      STEPSIZEFLAG = 0,
      SAMPLING_WINDOW = 1,
      EACHCLOCK = 500
    )
    
    set.seed(x$stoc_seed)
    suppressMessages(
      mod_obj <- fit_isingGraph3(
        DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
        METHOD = 'standard',
        CPP_CONTROL = cpp_ctrl,
        INIT = par_init,
        VERBOSEFLAG = 0
      )
    )
    
    
    mod_obj$lab <- lab
    mod_obj
  }), cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))



#### BERNOULLI ####
load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'bernoulli'
path <- paste0('isi/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
  purrr::transpose(.l = sim_settings ),
  function(x)try({
    
    n <- x$sample_size
    set.seed(x$stoc_seed)
    data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
    
    subtraj <- seq(0, x$maxiter, 250)
    cpp_ctrl <- list(
      MAXT = x$maxiter,
      BURN = x$burn,
      STEPSIZE = x$stepsize,
      SCALEVEC = x$scls,
      PAR1 = 1,
      PAR2 = x$stepsize,
      PAR3 = .501,
      NU = 1,
      STEPSIZEFLAG = 0,
      SAMPLING_WINDOW = 1,
      EACHCLOCK = 500
    )
    
    set.seed(x$stoc_seed)
    suppressMessages(
      mod_obj <- fit_isingGraph3(
        DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
        METHOD = 'bernoulli',
        CPP_CONTROL = cpp_ctrl,
        INIT = par_init,
        VERBOSEFLAG = 0
      )
    )
    
    mod_obj$lab <- lab
    mod_obj
  }), cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))


#### HYPER ####
load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'hyper'
path <- paste0('isi/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
  purrr::transpose(.l = sim_settings ),
  function(x)try({
    
    n <- x$sample_size
    set.seed(x$stoc_seed)
    data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
    
    subtraj <- seq(0, x$maxiter, 250)
    cpp_ctrl <- list(
      MAXT = x$maxiter,
      BURN = x$burn,
      STEPSIZE = x$stepsize,
      SCALEVEC = x$scls,
      PAR1 = 1,
      PAR2 = x$stepsize,
      PAR3 = .501,
      NU = 1,
      STEPSIZEFLAG = 0,
      SAMPLING_WINDOW = 1,
      EACHCLOCK = 500
    )
    
    set.seed(x$stoc_seed)
    suppressMessages(
      mod_obj <- fit_isingGraph3(
        DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
        METHOD = 'hyper',
        CPP_CONTROL = cpp_ctrl,
        INIT = par_init,
        VERBOSEFLAG = 0
      )
    )
    
    mod_obj$lab <- lab
    mod_obj
  }), cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))

#### RECYCLE STANDARD ####
load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'recycle_standard'
path <- paste0('isi/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
  purrr::transpose(.l = sim_settings ),
  function(x)try({
    
    n <- x$sample_size
    set.seed(x$stoc_seed)
    data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
    
    subtraj <- seq(0, x$maxiter, 250)
    cpp_ctrl <- list(
      MAXT = x$maxiter,
      BURN = x$burn,
      STEPSIZE = x$stepsize,
      SCALEVEC = x$scls,
      PAR1 = 1,
      PAR2 = x$stepsize,
      PAR3 = .501,
      NU = 1,
      STEPSIZEFLAG = 0,
      SAMPLING_WINDOW = 100,
      EACHCLOCK = 500
    )
    
    set.seed(x$stoc_seed)
    suppressMessages(
      mod_obj <- fit_isingGraph3(
        DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
        METHOD = 'recycle_standard',
        CPP_CONTROL = cpp_ctrl,
        INIT = par_init,
        VERBOSEFLAG = 0
      )
    )
    
    
    mod_obj$lab <- lab
    mod_obj
  }), cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))



#### RECYCLE HYPER ####
load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'recycle_hyper'
path <- paste0('isi/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')

est_obj_list <- pbapply::pblapply(
  purrr::transpose(.l = sim_settings ),
  function(x)try({
    
    n <- x$sample_size
    set.seed(x$stoc_seed)
    data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
    
    subtraj <- seq(0, x$maxiter, 250)
    cpp_ctrl <- list(
      MAXT = x$maxiter,
      BURN = x$burn,
      STEPSIZE = x$stepsize,
      SCALEVEC = x$scls,
      PAR1 = 1,
      PAR2 = x$stepsize,
      PAR3 = .501,
      NU = 1,
      STEPSIZEFLAG = 0,
      SAMPLING_WINDOW = 100,
      EACHCLOCK = 500
    )
    
    set.seed(x$stoc_seed)
    suppressMessages(
      mod_obj <- fit_isingGraph3(
        DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
        METHOD = 'recycle_hyper',
        CPP_CONTROL = cpp_ctrl,
        INIT = par_init,
        VERBOSEFLAG = 0
      )
    )
    
    mod_obj$lab <- lab
    mod_obj
  }), cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))

#### COLLECT #####
load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))

true_tib <- tibble(par = 1:length(true_theta), true_val = true_theta)
load(paste0('isi/data_files/', settingLab, '/ucminf/est.Rda'))
num_est <- sim_settings %>%
  mutate(stepsize = NA, scaling = NA, scls = NA, maxiter = NA, burn = NA, meth = 'ucminf') %>%
  distinct() %>%
  mutate(mod = est_obj_list)

load(paste0('isi/data_files/', settingLab, '/standard/est.Rda'))
standard_est <- sim_settings %>%
  mutate(meth = 'standard') %>%
  distinct() %>%
  mutate(mod = est_obj_list)

load(paste0('isi/data_files/', settingLab, '/bernoulli/est.Rda'))
bernoulli_est <- sim_settings %>%
  mutate(meth = 'bernoulli') %>%
  distinct() %>%
  mutate(mod = est_obj_list)

load(paste0('isi/data_files/', settingLab, '/hyper/est.Rda'))
hyper_est <- sim_settings %>%
  mutate(meth = 'hyper') %>%
  distinct() %>%
  mutate(mod = est_obj_list)

load(paste0('isi/data_files/', settingLab, '/recycle_standard/est.Rda'))
recycle_standard_est <- sim_settings %>%
  mutate(meth = 'recycle_standard') %>%
  distinct() %>%
  mutate(mod = est_obj_list)

load(paste0('isi/data_files/', settingLab, '/recycle_hyper/est.Rda'))
recycle_hyper_est <- sim_settings %>%
  mutate(meth = 'recycle_hyper') %>%
  distinct() %>%
  mutate(mod = est_obj_list)

stoc_est <- recycle_standard_est %>%
  bind_rows(bernoulli_est) %>%
  bind_rows(hyper_est) %>%
  bind_rows(standard_est) %>%
  bind_rows(recycle_hyper_est)

path_tab <- stoc_est %>%
  mutate(path_theta = map(mod, ~try(get_tidy_path3(.x, 'path_av_theta')))) %>%
  select(-mod)

path_tab |> filter(row_number()==46)
save(num_est, stoc_est, path_tab, true_tib, file = paste0('isi/data_files/', settingLab, '/est.Rda'))
######## variance ####

load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))
load(paste0('isi/data_files/',settingLab,'/est.Rda'))
selected_passes <- seq(.5,npasses,.5)
ncores <- 1

#### Computation ####
path_selected <- path_tab %>%
  unnest(path_theta) %>%
  mutate(pass = iter/sample_size) %>%
  filter(pass %in% selected_passes) %>%
  distinct()

cat('Variance stochastic estimators:\n')
var_list <- pbapply::pblapply(path_selected %>% purrr::transpose(),
                              function(x)try({
                                
                                # regenerate data
                                n <- x$sample_size
                                set.seed(x$stoc_seed)
                                data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
                                
                                # read method
                                method <- switch (x$meth,
                                                  'GD'               = 0,
                                                  'ucminf'           = 0,
                                                  'standard'         = 1,
                                                  'recycle_standard' = 1,
                                                  'bernoulli'        = 2,
                                                  'hyper'            = 2,
                                                  'recycle_hyper'    = 2)
                                
                                out <- sampleVar(
                                  THETA = x$path_av_theta,
                                  DATA = as.matrix(data),
                                  CONSTRAINTS = Q,
                                  NU = 1,
                                  METHOD = method,
                                  RANGE = x$iter - x$burn,
                                  TOTFLAG = T,
                                  PRINTFLAG = F
                                  
                                )
                                
                                se <- out$se
                                
                                
                                vt_tot <- out$var$var_tot %>% diag() %>% sum()
                                vt_stoc <- out$var$var_stoc %>% diag() %>% sum()
                                vt_stat <- out$var$var_stoc %>% diag() %>% sum()
                                tr <- tibble( var_type = c('stoc', 'stat', 'tot'), tr = c(vt_stoc, vt_stat, vt_tot))
                                
                                
                                return(list('se'=se, 'trace' = tr ))
                              }),
                              cl = ncores)

path_selected$asy_var <- var_list
save(path_selected, file = paste0('isi/data_files/', settingLab,'/','path_variance.Rda'))


cat('Variance numeric estimator:\n')
var_list_num <- pbapply::pblapply(num_est %>% purrr::transpose(),
                                  function(x)try({
                                    
                                    # regenerate data
                                    n <- x$sample_size
                                    set.seed(x$stoc_seed)
                                    data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
                                    
                                    # read method
                                    method <- switch (x$meth,
                                                      'GD'     = 0,
                                                      'ucminf' = 0)
                                    
                                    out <- sampleVar(
                                      THETA = x$mod$theta,
                                      DATA = as.matrix(data),
                                      CONSTRAINTS = Q,
                                      NU = 1,
                                      METHOD = method,
                                      RANGE = 1,
                                      TOTFLAG = T,
                                      PRINTFLAG = F
                                    )
                                    
                                    se <- out$se
                                    
                                    vt_tot <- out$var$var_tot %>% diag() %>% sum()
                                    vt_stoc <- out$var$var_stoc %>% diag() %>% sum()
                                    vt_stat <- out$var$var_stoc %>% diag() %>% sum()
                                    tr <- tibble( var_type = c('stoc', 'stat', 'tot'), tr = c(vt_stoc, vt_stat, vt_tot))
                                    
                                    return(list('se'=se, 'trace' = tr ))
                                  }), cl = ncores)
num_est$asy_var <- var_list_num

save(num_est, file = paste0('isi/data_files/', settingLab,'/','num_variance.Rda'))
