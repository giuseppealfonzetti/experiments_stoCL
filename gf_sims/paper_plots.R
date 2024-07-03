library(tidyverse)
library(ggh4x)
plots_path <- paste0('gf_sims/plots/')

colzA <- scales::viridis_pal(option = 'G')(8)
colzB <- scales::viridis_pal(option = 'B')(8)
colzC <- scales::viridis_pal(option = 'A')(8)
colz <- c(colzA[c(3, 4, 6)], colzB[c(6, 4, 7)])
colz1 <- colz[c(3,4,5)]
numcol <- '#0B0405FF'


#### Mean square error ####
mse_tot_20 <- readRDS("gf_sims/plots/I/gg/mse_tot.rds")
mse_tot_30 <- readRDS("gf_sims/plots/II/gg/mse_tot.rds")

join_stoc_mse <- mse_tot_20$data %>%
    mutate(p = 20) %>%
    bind_rows(
        mse_tot_30$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size) #

join_num_mse <- mse_tot_20$layers[[2]]$data %>%
    mutate(p = 20) %>%
    bind_rows(
        mse_tot_30$layers[[2]]$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size)

gg_join_mse <- join_stoc_mse %>%
    filter(meth %in% c('recycle_hyper_500', 'recycle_standard_500')) %>%
    # mutate(meth = map_chr(meth, ~if_else(.x=='recycle_hyper_500', 'recycle_hyper', .x))) %>%
    ggplot(aes(x = pass, y = log(mse)))+
    geom_rect(aes(xmin=-Inf,xmax=.25,ymin=-Inf,ymax=Inf),
              fill="lightgrey", color =NA)+
    geom_hline(data = join_num_mse %>%  filter(par_type == 'all') ,
               aes(yintercept = log(mse)), col = numcol, linetype = 'dashed', alpha = 0.5)+
    geom_line(aes(col = factor(meth)), alpha = .8, linewidth = 1)+
    facet_nested(p~n, scales = 'free', labeller = purrr::partial(label_both, sep = " = ")) +
    theme_bw()+
    labs(col = 'Sampling scheme:', x = 'Number of passes through the data', y = 'Log mean square error') +
    scale_color_manual(values = colz1) +
    theme(legend.position = 'bottom')
gg_join_mse
ggsave(plot = gg_join_mse, filename = paste0(plots_path, 'gg_join_mse.pdf'), width = 8, height = 5)


#### Variance trace ####
var_tr_20 <- readRDS("gf_sims/plots/I/gg/var_tr.rds")
var_tr_30 <- readRDS("gf_sims/plots/II/gg/var_tr.rds")

join_stoc_tr <- var_tr_20$data %>%
    mutate(p = 20) %>%
    bind_rows(
        var_tr_30$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size)

join_num_tr <- var_tr_20$layers[[1]]$data %>%
    mutate(p = 20) %>%
    bind_rows(
        var_tr_30$layers[[1]]$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size)

gg_join_var <- join_stoc_tr %>%
    filter(meth %in% c('recycle_hyper_500', 'recycle_standard_500')) %>%
    # mutate(meth = map_chr(meth, ~if_else(.x=='recycle_hyper_500', 'recycle_hyper', .x))) %>%
    ggplot(aes(x = pass, y = (tr)))+
    geom_hline(data = join_num_tr,
               aes(yintercept = (tr)),  col = numcol, linetype = 'dashed', alpha = 0.5) +
    geom_line(aes(col = meth), linewidth = 1, alpha = .8)+
    facet_nested(p~n, scales = 'free', labeller = purrr::partial(label_both, sep = " = ")) +
    theme_bw() +
    labs(x = 'Number of passes through the data', y = 'Variance trace', col = 'Sampling scheme:')+
    scale_color_manual(values = colz1) +
    theme(legend.position = 'bottom')

ggsave(plot = gg_join_var, filename = paste0(plots_path, 'gg_join_var.pdf'), width = 8, height = 5)

#### Inference ####
cov_lines20 <- readRDS("gf_sims/plots/I/gg/cov_lines.rds")
cov_lines30 <- readRDS("gf_sims/plots/II/gg/cov_lines.rds")

join_stoc_cov <- cov_lines20$data %>%
    mutate(p = 20) %>%
    bind_rows(
        cov_lines30$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size)

join_num_cov <- cov_lines20$layers[[3]]$data %>%
    mutate(p = 20) %>%
    bind_rows(
        cov_lines30$layers[[3]]$data %>%
            mutate(p = 30)
    ) %>%
    mutate(npasses = 3.05) %>%
    rename(n = sample_size)

gg_join_cov <- join_stoc_cov %>%
    filter(meth %in% c('recycle_hyper_500', 'recycle_standard_500')) %>%
    # mutate(meth = map_chr(meth, ~if_else(.x=='recycle_hyper_500', 'recycle_hyper', .x))) %>%
    mutate(n = factor(n, levels = unique(join_num_cov$n), labels = paste0('n = ', unique(join_num_cov$n))),
           p = factor(p, levels = unique(join_num_cov$p), labels = paste0('p = ', unique(join_num_cov$p)))) %>%
    ggplot(aes(x = pass, y = coverage, col = meth, group = interaction(meth, par)) )+
    geom_hline(yintercept = c(.9,.95, 1), linetype = 'dashed', alpha = .5) +
    geom_line(alpha = .8)+
    geom_point(data = join_num_cov  %>%
                   mutate(n = factor(n, levels = unique(join_num_cov$n), labels = paste0('n = ', unique(join_num_cov$n))),
                          p = factor(p, levels = unique(join_num_cov$p), labels = paste0('p = ', unique(join_num_cov$p)))) , alpha = .5, col = numcol, shape = 4, size = 2) +
    facet_nested( p ~ n + regime, scales = 'free') +
    theme_bw() +
    labs(y = 'Coverage', x = 'Number of passes through the data', col = 'Sampling scheme:') +
    scale_color_manual(values = colz1) +
    theme(legend.position = 'bottom')
gg_join_cov
ggsave(plot = gg_join_cov, filename = paste0(plots_path, 'gg_join_cov.pdf'), width = 8, height = 5)


#### Appendix plots ####
appendix_mse_20 <- readRDS("gf_sims/plots/I/gg/appendix_mse.rds")
appendix_mse_30 <- readRDS("gf_sims/plots/II/gg/appendix_mse.rds")

join_appendix_stoc_mse <- appendix_mse_20$data %>%
    mutate(p = 20) %>%
    bind_rows(
        appendix_mse_30$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size) %>%
    mutate(meth = map_chr(meth, ~if_else(.x=='recycle_standard_500', 'recycle_standard', .x)))

join_appendix_num_mse <- appendix_mse_20$layers[[2]]$data %>%
    mutate(p = 20) %>%
    bind_rows(
        appendix_mse_30$layers[[2]]$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size)

join_appendix_stoc_mse %>%
    ggplot(aes(x = pass, y = log(mse)))+
    geom_rect(aes(xmin=-Inf,xmax=.25,ymin=-Inf,ymax=Inf),
              fill="lightgrey", color =NA)+
    geom_hline(data = join_num_mse %>% ungroup() %>% select(-stepsize) %>%  expand_grid(stepsize = unique(join_appendix_stoc_mse$stepsize)),
               aes(yintercept = log(mse)), col = numcol, linetype = 'dashed', alpha = 0.5)+
    geom_line(aes(col = factor(meth)), alpha = .8) +
    facet_nested(par_type+p~n+stepsize, scales = 'free') +
    theme_bw()+
    labs(col = 'Sampling scheme:', x = 'Number of passes through the data', y = 'Log mean square error') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')

join_appendix_stoc_mse %>%
    filter(meth == 'recycle_standard', n == 10000, par_type == 'all', p == 20, pass == 3)

gg_app_join_stepsize <- join_appendix_stoc_mse %>%
    filter(par_type == 'all', pass %in% c( 1,2, 3)) %>%
    ggplot(aes( x = stepsize, y = log(mse), col = meth)) +
    geom_line(alpha = .8)+
    geom_point(alpha =.5) +
    facet_nested(n~p+pass, labeller = purrr::partial(label_both, sep = " = "), scales = 'free') +
    theme_bw() +
    labs(col = 'Sampling scheme:', x = 'Stepsize', y = 'Log mean square error') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')
ggsave(plot = gg_app_join_stepsize, filename = paste0(plots_path, 'gg_app_join_stepsize.pdf'), width = 10, height = 8)

gg_app_join_mse <- join_appendix_stoc_mse %>%
    mutate(n = factor(n, levels = unique(join_appendix_stoc_mse$n), labels = paste0('n = ', unique(join_appendix_stoc_mse$n))),
           p = factor(p, levels = unique(join_appendix_stoc_mse$p), labels = paste0('p = ', unique(join_appendix_stoc_mse$p)))) %>%
    filter(par_type!='all') %>%
    ggplot(aes(x = pass, y = log(mse)))+
    geom_rect(aes(xmin=-Inf,xmax=.25,ymin=-Inf,ymax=Inf),
              fill="lightgrey", color =NA)+
    geom_hline(data = join_appendix_num_mse %>% filter(par_type != 'all') %>%  ungroup() %>% select(-stepsize) %>%  expand_grid(stepsize = unique(join_appendix_stoc_mse$stepsize)) %>%
                   mutate(n = factor(n, levels = unique(join_appendix_stoc_mse$n), labels = paste0('n = ', unique(join_appendix_stoc_mse$n))),
                          p = factor(p, levels = unique(join_appendix_stoc_mse$p), labels = paste0('p = ', unique(join_appendix_stoc_mse$p)))),
               aes(yintercept = log(mse)), col = numcol, linetype = 'dashed', alpha = 0.5)+
    geom_line(aes(col = factor(meth)), alpha = .8) +
    facet_nested(n+stepsize ~ p+par_type, scales = 'free') +
    theme_bw()+
    labs(col = 'Sampling scheme:', x = 'Number of passes through the data', y = 'Log mean square error') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')

ggsave(plot = gg_app_join_mse, filename = paste0(plots_path, 'gg_app_join_mse.pdf'), width = 10, height = 10)

gg_app_zoom <- join_appendix_stoc_mse %>%
    filter(pass > 0.25, meth != 'recycle_standard', stepsize == 2) %>%
    mutate(n = factor(n, levels = unique(join_appendix_stoc_mse$n), labels = paste0('n = ', unique(join_appendix_stoc_mse$n))),
           p = factor(p, levels = unique(join_appendix_stoc_mse$p), labels = paste0('p = ', unique(join_appendix_stoc_mse$p)))) %>%
    filter(par_type=='all') %>%
    ggplot(aes(x = pass, y = log(mse)))+
    geom_line(aes(col = factor(meth)), alpha = .8) +
    facet_nested(n ~ p, scales = 'free') +
    theme_bw()+
    labs(col = 'Sampling scheme:', x = 'Number of passes through the data', y = 'Log mean square error') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')
ggsave(plot = gg_app_zoom, filename = paste0(plots_path, 'gg_app_zoom.pdf'), width = 10, height = 5)


#### variance ####
appendix_var_tr_20 <- readRDS("gf_sims/plots/I/gg/appendix_var_tr.rds")
appendix_var_tr_30 <- readRDS("gf_sims/plots/II/gg/appendix_var_tr.rds")

join_stoc_tr <- appendix_var_tr_20$data %>%
    mutate(p = 20) %>%
    bind_rows(
        appendix_var_tr_30$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size, step = stepsize)

join_num_tr <- appendix_var_tr_20$layers[[1]]$data %>%
    mutate(p = 20) %>%
    bind_rows(
        appendix_var_tr_30$layers[[1]]$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size, step = stepsize)

join_stoc_tr  %>%
    ggplot(aes(x = pass, y = (tr)))+
    geom_hline(data = join_num_tr,
               aes(yintercept = (tr)),  col = numcol, linetype = 'dashed', alpha = 0.5) +
    geom_line(aes(col = meth), linewidth = 1, alpha = .8)+
    facet_nested(p+step~n, scales = 'free', labeller = purrr::partial(label_both, sep = " = ")) +
    theme_bw() +
    labs(x = 'Number of passes through the data', y = 'Variance trace', col = 'Sampling scheme:')+
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')

#### coverage ####
appendix_cov_lines_20 <- readRDS("gf_sims/plots/I/gg/appendix_cov_lines.rds")
appendix_cov_lines_30 <- readRDS("gf_sims/plots/II/gg/appendix_cov_lines.rds")

join_stoc_cov <- appendix_cov_lines_20$data %>%
    mutate(p = 20) %>%
    bind_rows(
        appendix_cov_lines_30$data %>%
            mutate(p = 30)
    ) %>%
    rename(n = sample_size, step = stepsize)

join_num_cov <- appendix_cov_lines_20$layers[[3]]$data %>%
    mutate(p = 20) %>%
    bind_rows(
        appendix_cov_lines_30$layers[[3]]$data %>%
            mutate(p = 30)
    ) %>%
    mutate(npasses = 3.05) %>%
    rename(n = sample_size)

gg_app_join_coverage <- join_stoc_cov  %>%
    mutate(n = factor(n, levels = unique(join_num_cov$n), labels = paste0('n = ', unique(join_num_cov$n))),
           p = factor(p, levels = unique(join_num_cov$p), labels = paste0('p = ', unique(join_num_cov$p)))) %>%
    ggplot(aes(x = pass, y = coverage, col = meth, group = interaction(meth, par)) )+
    geom_hline(yintercept = c(.9,.95, 1), linetype = 'dashed', alpha = .5) +
    geom_line(alpha = .8)+
    geom_point(data = join_num_cov %>%  expand_grid(step = unique(join_stoc_cov$step))  %>%
                   mutate(n = factor(n, levels = unique(join_num_cov$n), labels = paste0('n = ', unique(join_num_cov$n))),
                          p = factor(p, levels = unique(join_num_cov$p), labels = paste0('p = ', unique(join_num_cov$p)))) , alpha = .5, col = numcol, shape = 4, size = 2) +
    facet_nested( p + step ~ n + regime, scales = 'free') +
    theme_bw() +
    labs(y = 'Coverage', x = 'Number of passes through the data', col = 'Sampling scheme:') +
    scale_color_manual(values = colz) +
    theme(legend.position = 'bottom')

ggsave(plot = gg_app_join_coverage, filename = paste0(plots_path, 'gg_app_join_coverage.pdf'), width = 10, height = 10)

