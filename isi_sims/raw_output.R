### dependencies ####
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)
library(tidyverse)

#### choose setting folder ####
settingLab <- 'I'
load(paste0('isi/data_files/',settingLab,'/setting_init.Rda'))
load(paste0('isi/data_files/',settingLab,'/sim_setup.Rda'))
plots_path <- paste0('isi/data_files/', settingLab, '/plots/')
if(!dir.exists(plots_path)){dir.create(plots_path); cat('[Create setting-related data folder]')}
if(!dir.exists(paste0(plots_path,'gg/'))){dir.create(paste0(plots_path,'gg/')); cat('[Create setting-related data folder]')}

colzA <- scales::viridis_pal(option = 'G')(8)
colzB <- scales::viridis_pal(option = 'B')(8)
colzC <- scales::viridis_pal(option = 'A')(8)

barplot(rep(1, length(colzA)), col = colzA)
barplot(rep(1, length(colzB)), col = colzB)
barplot(rep(1, length(colzC)), col = colzC)
colz <- c(colzA[c(3, 4, 6)], colzB[c(6, 4, 7)])
#colz <- colz[c(3,4,5)]

barplot(rep(1, length(colz)), col = colz)
numcol <- '#0B0405FF'

#### pointwise results ####
load(paste0('isi/data_files/',settingLab,'/est.Rda'))
name_par <- function(par){
    if(par > p){'edge coefficients'}
    else{'intercepts'}
}


# stochastic parameters trajectories: mse
par_path <- path_tab %>%
    unnest(path_theta) %>%
    mutate(intercepts = map_dbl(path_av_theta, ~mean((.x[1:p]-true_theta[1:p])^2)),
           edges = map_dbl(path_av_theta, ~mean((.x[(p+1):d]-true_theta[(p+1):d])^2)),
           all = map_dbl(path_av_theta, ~mean((.x-true_theta)^2))) %>%
    gather(key = 'par_type', value = 'mse', intercepts, edges, all) %>%
    group_by(sample_size, meth, stepsize, scaling, iter, par_type) %>%
    summarise(mse = mean(mse))



# numerical estimates as comparison: mse
num_par <- num_est %>%
    mutate(theta = map(mod, ~.x$theta)) %>%
    select(-mod) %>%
    mutate(intercepts = map_dbl(theta, ~mean((.x[1:p]-true_theta[1:p])^2)),
           edges = map_dbl(theta, ~mean((.x[(p+1):d]-true_theta[(p+1):d])^2)),
           all = map_dbl(theta, ~mean((.x-true_theta)^2))) %>%
    gather(key = 'par_type', value = 'mse', intercepts, edges, all) %>%
    group_by(sample_size, meth, stepsize, scaling, par_type) %>%
    summarise(mse = mean(mse))

appendix_mse <- par_path %>%
    mutate(pass = iter/sample_size) %>%
    ggplot(aes(x = pass, y = log(mse)))+
    geom_rect(aes(xmin=-Inf,xmax=burn,ymin=-Inf,ymax=Inf),
              fill="lightgrey", color =NA)+
    geom_hline(data = num_par %>%
                   ungroup() %>%
                   select(-stepsize, -meth) %>%
                   expand_grid(stepsize=unique(stoc_est$stepsize), meth = unique(stoc_est$meth)),
               aes(yintercept = log(mse)), col = numcol, linetype = 'dashed', alpha = 0.5) +
    geom_line(aes(col = factor(meth)), alpha =.8) +
    facet_nested(par_type~sample_size+stepsize, scales = 'free') +
    theme_bw()+
    labs(col = 'Sampling scheme:', x = 'Number of passes through the data', y = 'Log mean square error') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')
appendix_mse
saveRDS(appendix_mse, file = paste0(plots_path,'gg/appendix_mse.rds'))
ggsave(plot = appendix_mse, filename = paste0(plots_path, 'appendix_mse.pdf'), width = 10, height = 5)

mse <- par_path %>%
    filter(stepsize == 1) %>%
    mutate(pass = iter/sample_size) %>%
    ggplot(aes(x = pass, y = log(mse)))+
    geom_rect(aes(xmin=-Inf,xmax=burn,ymin=-Inf,ymax=Inf),
              fill="lightgrey", color =NA)+
    geom_hline(data = num_par %>%
                   #filter(par_type != 'all') %>%
                   ungroup() %>%
                   select(-stepsize, -meth) %>%
                   expand_grid( meth = unique(stoc_est$meth)),
               aes(yintercept = log(mse)), col = numcol, linetype = 'dashed', alpha = 0.5) +
    geom_line(aes(col = factor(meth)), alpha = .8, linewidth = 1) +
    facet_nested(par_type~sample_size, scales = 'free') +
    theme_bw()+
    labs(col = 'Sampling scheme:', x = 'Number of passes through the data', y = 'Log mean square error') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')
mse
saveRDS(mse, file = paste0(plots_path,'gg/mse.rds'))
ggsave(plot = mse, filename = paste0(plots_path, 'mse.pdf'), width = 8, height = 4)

mse_tot <- par_path %>%
    filter(stepsize == 1, par_type == 'all') %>%
    mutate(pass = iter/sample_size) %>%
    ggplot(aes(x = pass, y = log(mse)))+
    geom_rect(aes(xmin=-Inf,xmax=burn,ymin=-Inf,ymax=Inf),
              fill="lightgrey", color =NA)+
    geom_hline(data = num_par %>%
                   filter(par_type == 'all') %>%
                   ungroup() %>%
                   select(-stepsize, -meth) %>%
                   expand_grid(stepsize=unique(stoc_est$stepsize), meth = unique(stoc_est$meth)),
               aes(yintercept = log(mse)), col = numcol, linetype = 'dashed', alpha = 0.5) +

    geom_line(aes(col = factor(meth)), alpha = .8, linewidth = 1) +
    facet_nested(~sample_size, scales = 'free') +
    theme_bw()+
    labs(col = 'Sampling scheme:', x = 'Number of passes through the data', y = 'Log mean square error') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')
mse_tot
saveRDS(mse_tot, file = paste0(plots_path,'gg/mse_tot.rds'))
ggsave(plot = mse_tot, filename = paste0(plots_path, 'mse_tot.pdf'), width = 8, height = 4)

#### Observed variance trace ####
selected_passes <- seq(.5,npasses,.5)

# tajectory stochastic parameters: average across simulations
sim_average <- path_tab %>%
    unnest(path_theta) %>%
    mutate(pass = iter/sample_size) %>%
    filter(pass %in% selected_passes) %>%
    select(-scls, - iter) %>%
    group_by(sample_size, meth, stepsize, pass, scaling) %>%
    summarise(av = list(colSums(purrr::reduce(path_av_theta, rbind))/n()))

# tajectory stochastic parameters: variance trace across simulations
sim_trace <- path_tab %>%
    unnest(path_theta) %>%
    mutate(pass = iter/sample_size) %>%
    filter(pass %in% selected_passes) %>%
    left_join(sim_average, by = c('sample_size', 'meth', 'stepsize', 'pass', 'scaling')) %>%
    mutate(sq = map2_dbl(path_av_theta, av, ~sum((.x-.y)^2))) %>%
    select(-scls, -path_av_theta, -av, -iter, -maxiter) %>%
    group_by(sample_size, meth, stepsize, pass, scaling) %>%
    summarise(tr = sum(sq)/(n()-1))

# numerical estimates: average across simulations
num_sim_average <- num_est %>%
    mutate(theta = map(mod, ~.x$theta)) %>%
    select(-mod) %>%
    group_by(sample_size, meth, stepsize) %>%
    summarise(av = list(colSums(purrr::reduce(theta, rbind))/n()))

# # numerical estimates: variance trace across simulations
num_sim_variance <- num_est %>%
    mutate(theta = map(mod, ~.x$theta)) %>%
    select(-mod) %>%
    left_join(num_sim_average, by = c('sample_size', 'meth', 'stepsize')) %>%
    mutate(sq = map2_dbl(theta, av, ~sum((.x-.y)^2))) %>%
    select(-scls, -theta, -av, -maxiter) %>%
    group_by(sample_size, meth) %>%
    summarise(tr = sum(sq)/(n()-1))

appendix_var_tr <- sim_trace %>%
    ggplot(aes(x = pass, y = (tr)))+
    geom_hline(data = num_sim_variance %>%
                   expand_grid(stepsize=unique(stoc_est$stepsize)),
               aes(yintercept = (tr)), col = numcol, linetype = 'dashed', alpha = .5) +
    geom_line(aes(col = meth), linewidth = 1, alpha = .6)+
    facet_nested(stepsize~sample_size, scales = 'free') +
    theme_bw() +
    labs(x = 'Number of passes through the data', y = 'Variance trace', col = 'Sampling scheme:')+
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')

saveRDS(appendix_var_tr, file = paste0(plots_path,'gg/appendix_var_tr.rds'))
ggsave(plot = appendix_var_tr, filename = paste0(plots_path, 'appendix_var_tr.pdf'), width = 10, height = 5)

var_tr <- sim_trace %>%
    filter(stepsize == 1) %>%
    ggplot(aes(x = pass, y = (tr)))+
    geom_hline(data = num_sim_variance %>%
                   expand_grid(stepsize=unique(stoc_est$stepsize)),
               aes(yintercept = (tr)), col = numcol, linetype = 'dashed', alpha = .5) +
    geom_line(aes(col = meth), linewidth = 1, alpha = .6)+
    facet_nested(~sample_size, scales = 'free') +
    theme_bw() +
    labs(x = 'Number of passes through the data', y = 'Variance trace', col = 'Sampling scheme:')+
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')
var_tr
saveRDS(var_tr, file = paste0(plots_path,'gg/var_tr.rds'))
ggsave(plot = var_tr, filename = paste0(plots_path, 'var_tr.pdf'), width = 8, height = 3)

#### time ####
clock_tib <- num_est %>%
    mutate(
        clock = 'main',
        time = map_dbl(mod,~.x$time)) %>%
    select(-mod, -scls) %>%
    bind_rows(
        stoc_est %>%
            mutate(
                main = map_dbl(mod,~.x$clock[1,2]),
                sampling = map_dbl(mod,~.x$clock[2,2]),
                gradient = map_dbl(mod,~.x$clock[3,2]),
                update = map_dbl(mod,~.x$clock[4,2])
            ) %>%
            select(-mod, -scls) %>%
            gather(key = 'clock', value = 'time', main, sampling, gradient, update)
    ) %>%
    mutate(clock = factor(clock, levels = c('main', 'sampling', 'gradient', 'update'), ordered = T))
gg_clock_st <- clock_tib %>%
    filter(clock!='main') %>%
    ggplot(aes(x = meth, y = (time), fill = meth)) +
    geom_boxplot() +
    facet_grid(clock~sample_size, scales = 'free') +
    theme_bw()+
    labs(x = ' ')+
    theme(axis.text.x = element_text(angle = 90), legend.position = 'bottom' )+
    scale_fill_viridis_d()

clock_tib %>%
    group_by(sample_size, clock, meth) %>%
    summarise(q25=quantile(time, .25),
              q50=quantile(time, .5),
              q75=quantile(time, .75)) %>%
    print(n=100)

clock_table <- clock_tib %>%
    group_by(sample_size, clock, meth) %>%
    summarise(time = mean(time)) %>%
    spread(meth, time) %>%
    mutate(clock = factor(clock, levels = c('main', 'sampling', 'gradient', 'update'),
                          labels = c('Total', 'Sampling Step', 'Approximation Step', 'Update Step')))%>%
    rename(n = sample_size, operation = clock) %>%
    print(n=100)

library(xtable)
print(xtable(clock_table,
             display = rep('e', ncol(clock_table) + 1)),
      include.rownames=FALSE,
      file = paste0(plots_path, 'clock.tex'),
      math.style.exponents = TRUE,
      booktabs = T,
      size = 'footnotesize')


gg_clock_main <- clock_tib %>%
    filter(clock=='main') %>%
    ggplot(aes(x = meth, y = (time), fill = meth)) +
    geom_boxplot(aes(fill = NULL)) +
    facet_grid(clock~sample_size, scales = 'free') +
    theme_bw()+
    labs(x = ' ')+
    theme(axis.text.x = element_text(angle = 90), legend.position = 'bottom' )
gg_clock_main
gg_clock <- ggpubr::ggarrange(gg_clock_main, gg_clock_st, nrow = 2, heights = c(1,2))
gg_clock
ggsave(plot = gg_clock, filename = paste0(plots_path, 'gg_clock.pdf'), width = 15, height = 15)



#### Inference results ####
load(paste0('isi/data_files/', settingLab,'/','path_variance.Rda'))
load(paste0('isi/data_files/', settingLab,'/','num_variance.Rda'))
selected_passes <- seq(.5,npasses,.5)
num_est %>%
    select(stoc_seed, sample_size, meth, asy_var) %>%
    mutate(
        sd_tot = lapply(asy_var, function(x){
            tib <- tibble(
                par = 1:length(x$se$se_tot),
                sd_tot = x$se$se_tot
            )
            tib
        })
    ) %>% pluck('asy_var',1)
sample_sd_num <- num_est %>%
    select(stoc_seed, sample_size, meth, asy_var) %>%
    mutate(
        sd_tot = lapply(asy_var, function(x){
            tib <- tibble(
                par = 1:length(x$se$se_tot),
                sd_tot = x$se$se_tot
            )
            tib
        })
    ) %>%
    select(-asy_var) %>%
    unnest(c(sd_tot))

num_coverage <- num_est %>%
    mutate(
        coverage_tot  = map2(mod, asy_var, ~ true_theta > .x$theta-1.96*.y$se$se_tot  & true_theta < .x$theta+1.96*.y$se$se_tot )
    ) %>%
    mutate(
        coverage_tot = lapply(coverage_tot, function(x){
            tib <- tibble(
                par = 1:length(x),
                cov_tot = x
            )
            tib
        })

    ) %>%
    unnest(c(coverage_tot))

st_coverage <- path_selected %>%
    mutate(
        coverage_tot  = map2(path_av_theta, asy_var, ~ true_theta > .x-1.96*.y$se$se_tot  & true_theta < .x+1.96*.y$se$se_tot ),
        coverage_stat = map2(path_av_theta, asy_var, ~ true_theta > .x-1.96*.y$se$se_stat & true_theta < .x+1.96*.y$se$se_stat),
        coverage_stoc = map2(path_av_theta, asy_var, ~ true_theta > .x-1.96*.y$se$se_stoc & true_theta < .x+1.96*.y$se$se_stoc)
    ) %>%
    select(-asy_var) %>%
    mutate(
        coverage_tot = lapply(coverage_tot, function(x){
            tib <- tibble(
                par = 1:length(x),
                cov_tot = x
            )
            tib
        }),
        coverage_stat = lapply(coverage_stat, function(x){
            tib <- tibble(
                cov_stat = x
            )
            tib
        }),
        coverage_stoc = lapply(coverage_stoc, function(x){
            tib <- tibble(
                cov_stoc = x
            )
            tib
        })
    ) %>%
    unnest(c(coverage_tot, coverage_stat, coverage_stoc))


sample_st_coverage <- st_coverage %>%
    select(sample_size, meth, stepsize, scaling, pass, par, cov_stat, cov_tot, cov_stoc) %>%
    group_by(sample_size, meth, stepsize, scaling, pass, par) %>%
    summarise_at( vars(starts_with('cov')), mean) %>%
    gather(key = 'regime', val = 'coverage',starts_with('cov')) %>%
    mutate(regime = factor(regime, levels = c('cov_stat', 'cov_stoc', 'cov_tot'), labels = c('Regime 1', 'Regime 2', 'Regime 3')))



sample_num_coverage <- num_coverage %>%
    select(sample_size, meth, par, cov_tot) %>%
    group_by(sample_size, meth,  par) %>%
    summarise_at( vars(starts_with('cov')), mean) %>%
    mutate(regime = 'Regime 3', pass = npasses+.1) %>%
    rename(coverage = cov_tot)

appendix_cov_lines <- sample_st_coverage %>%
    ggplot(aes(x= pass, y = coverage, col = meth, group = interaction(meth, par)) )+
    geom_hline(yintercept = c(.9,.95, 1), linetype = 'dashed', alpha = .5) +
    geom_line()+
    geom_point(data = sample_num_coverage, alpha = .5, col = numcol, shape = 4, size = 2) +
    facet_nested(stepsize ~ sample_size + regime) +
    theme_bw() +
    labs(y = 'Coverage', x = 'Number of passes through the data', col = 'Sampling scheme:') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')
appendix_cov_lines
saveRDS(appendix_cov_lines, file = paste0(plots_path,'gg/appendix_cov_lines.rds'))
ggsave(plot = appendix_cov_lines, filename = paste0(plots_path, 'appendix_cov_lines.pdf'), width = 8, height = 8)

cov_lines <- sample_st_coverage %>%
    filter(stepsize == 1) %>%
    ggplot(aes(x= pass, y = coverage, col = meth, group = interaction(meth, par)) )+
    geom_hline(yintercept = c(.9,.95, 1), linetype = 'dashed', alpha = .5) +
    geom_line(alpha = .8)+
    geom_point(data = sample_num_coverage, alpha = .5, col = numcol, shape = 4, size = 2) +
    facet_nested( ~ sample_size + regime) +
    theme_bw() +
    labs(y = 'Coverage', x = 'Number of passes through the data', col = 'Sampling scheme:') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')
cov_lines
plotly::ggplotly(cov_lines)
saveRDS(cov_lines, file = paste0(plots_path,'gg/cov_lines.rds'))
ggsave(plot = cov_lines, filename = paste0(plots_path, 'cov_lines.pdf'), width = 8, height = 4)

