#### dependencies ####
library(gammaFrailty)
library(RcppClock)
library(dplyr)
library(tidyr)
library(purrr)
library(pbapply)
pbo = pboptions(type="txt") 

settingLab <- 'II'
### Setting setup ####

load(paste0('gf_sims/data_files/',settingLab,'/setting_init.Rda'))

### Simulation setup ####

load(paste0('gf_sims/data_files/',settingLab,'/sim_setup.Rda'))

#'
#' ### Point estimation
#'

ncores <- 1

#### randomized ####
load(paste0('gf_sims/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('gf_sims/data_files/',settingLab,'/sim_setup.Rda'))

lab <- 'randomized'
path <- paste0('gf_sims/data_files/', settingLab,'/',lab)
if(!dir.exists(path)){dir.create(path); cat('[Create setting-method data folder]')}
cat(lab,':\n')
probs <- c(.5, .25, .1)

est_obj_list <- pbapply::pblapply(
  purrr::transpose(.l = sim_settings |> select(-c(stepsize, scls, scaling, maxiter, burn)) |> distinct() |> expand_grid(probs) ),
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
    
    suppressMessages(
      mod_obj <- fit_gammaFrailty2(
        DATA_LIST = list('DATA' = data, 'X' = X),
        METHOD = 'randomized',
        CPP_CONTROL = list(PROB = x$probs, SEED=x$stoc_seed),
        VERBOSEFLAG= 0,
        PAIRS_RANGE = p,
        INIT = par_init,
        ITERATIONS_SUBSET = NULL,
        STRUCT = 'COMPOUND'
      )
    )
    
    
    mod_obj$lab <- lab
    mod_obj
  }, cl = ncores)
save(est_obj_list, file = paste0(path,'/est.Rda'))

randomized_est <- sim_settings |>
  mutate(stepsize = NA, scaling = NA, scls = NA, maxiter = NA, burn = NA, meth = 'randomized') |>
  distinct() |>
  expand_grid(probs) |>
  mutate(mod = est_obj_list)

rand_par <- randomized_est |>
  mutate(theta = map(mod, ~.x$theta)) |>
  select(-mod) |>
  mutate(repar_xi = map_dbl(theta, ~mean((.x[1]-repar_theta[1])^2)),
         repar_rho = map_dbl(theta, ~mean((.x[2]-repar_theta[2])^2)),
         intercepts = map_dbl(theta, ~mean((.x[(3+m):d]-repar_theta[(3+m):d])^2)),
         all = map_dbl(theta, ~mean((.x-repar_theta)^2))) |>
  gather(key = 'par_type', value = 'mse',repar_xi, repar_rho, intercepts, all) |>
  group_by(n=sample_size, meth, stepsize, scaling, probs, par_type) |>
  summarise(mse = mean(mse)) |>
  filter(par_type=='all') |>
  mutate(p=p)

#### plots ####
library(tidyverse)
library(ggh4x)
library(ggrepel)
plots_path <- paste0('gf_sims/plots/')

colzA <- scales::viridis_pal(option = 'G')(8)
colzB <- scales::viridis_pal(option = 'B')(8)
colzC <- scales::viridis_pal(option = 'A')(8)
colz <- c(colzA[c(3, 4, 6)], colzB[c(6, 4, 7)])

colz1 <- colz[c(3,4,5)]
numcol <- '#0B0405FF'


#### Mean square error ####
mse_tot_30 <- readRDS(paste0('gf_sims/plots/',settingLab,'/gg/mse_tot.rds'))

join_stoc_mse <- mse_tot_30$data |>
  mutate(p = 30)|>
  rename(n = sample_size) #

join_num_mse <- mse_tot_30$layers[[2]]$data |>
  mutate(p = 30)|>
  rename(n = sample_size)

gg_join_mse <- join_stoc_mse |>
  filter(meth %in% c('recycle_hyper_500', 'recycle_standard')) |>
  mutate(meth = map_chr(meth, ~if_else(.x=='recycle_hyper_500', 'recycle_hyper', .x))) |>
  ggplot(aes(x = pass, y = log(mse)))+
  geom_rect(aes(xmin=-Inf,xmax=.25,ymin=-Inf,ymax=Inf),
            fill="lightgrey", color =NA)+
  geom_hline(data = join_num_mse |>  filter(par_type == 'all') ,
             aes(yintercept = log(mse)), col = numcol, linetype = 'dashed', alpha = 0.5)+
  geom_line(aes(col = factor(meth)), alpha = .8, linewidth = 1)+
  facet_nested(p~n, scales = 'free', labeller = purrr::partial(label_both, sep = " = ")) +
  theme_bw()+
  labs(col = 'Sampling scheme:', x = 'Number of passes through the data', y = 'Log mean square error') +
  scale_color_manual(values = colz1) +
  theme(legend.position = 'bottom')
gg <- gg_join_mse +
  geom_hline(data = rand_par |> filter(probs!=.75), aes(yintercept = log(mse)), color='darkgrey', linetype = 'dashed') +
  geom_label_repel(data = rand_par |> filter(probs!=.75), aes(x=3, y = log(mse), label = probs), color='darkgrey', nudge_y = 1, nudge_x = -1)

ggsave(plot = gg, filename = paste0(plots_path, 'gg_randomizedpl.pdf'), width = 10, height = 5)

randomized_est |>
  mutate(
    clock = 'main',
    time = map_dbl(mod,~.x$time)) %>%
  select(-mod, -scls) |>
  group_by(sample_size, clock, meth, probs) %>%
  summarise(q25=quantile(time, .25),
            q50=quantile(time, .5),
            q75=quantile(time, .75)) %>%
  print(n=100)
