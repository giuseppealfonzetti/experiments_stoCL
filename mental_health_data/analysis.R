devtools::install_github(
    "giuseppealfonzetti/stIsing@18cda99",
    ref="main",
    force = F
)
library(tidyverse)
library(ggplot2)
library(stIsing)
library(tictoc)

devtools::install_github("jaredhuling/jcolors")
                         
ising_graph_lik <- function(par){
    ncl(as.matrix(h.mat), par, rep(T, d), VERBOSEFLAG = F)$nll/nrow(h.mat)
}
ising_graph_gradient <- function(par){
    ncl(as.matrix(h.mat), par, rep(T, d), VERBOSEFLAG = F)$ngradient
}

#### Load synthetic dataset ####
# The sythetic dataset is generated from numerical estimates computed on the
# real dataset.
# Differently from what happened on real data, the negative log-likelihood
# decrease on the holdout partition stops the algorithm a bit earlier.
# However, results are generally very similar.
load('mental_health_data/data/mock_dataset.rda')

library(xtable)
print(xtable(dictionary %>% select(Item = label, Area = area, Description = txt) %>% as.data.frame(),
             booktabs = T,
             size = 'footnotesize',
             align= 'lllp{3in}'),         include.rownames=FALSE,
      file = 'dictionary.tex',
      tabular.environment="longtable", floating = F)

#### split train and holdout ####
seed <- 123
set.seed(seed)
synth_data$id <- 1:nrow(synth_data)
train <- synth_data %>% dplyr::sample_frac(.90)
holdout  <- dplyr::anti_join(synth_data, train, by = 'id')

##### ESTIMATION #####

# convert data to matrix format
mat <- as.matrix(train %>% select(-id))
h.mat <- as.matrix(holdout %>% select(-id))


n <- nrow(mat)        # obsevations
p <- ncol(mat)        # nodes
d <- p + p*(p-1)/2    # parameters

par_init <- rep(0, d) # initialisation

# estimate
{
    cpp_ctrl_init <- list(
        MAXT = n,
        BURN = .25*n,
        STEPSIZE = 10,
        SCALEVEC = rep(1, d),
        PAR1 = 1,
        PAR2 = 1,
        PAR3 = .501,
        NU = 1,
        STEPSIZEFLAG = 0,
        SAMPLING_WINDOW = n,
        SEED = 1,
        EACH = round(.25*n, 0),
        EACHCLOCK = n
    )
    tictoc::tic()
    tn <- stepsize_tuning4(
        DATA_LIST = list(DATA = mat, CONSTRAINTS = rep(T, d), HOLDOUT = h.mat),
        METHOD = 'recycle_hyper',
        CPP_CONTROL = cpp_ctrl_init,
        STEPSIZE_INIT = 5,
        LENGTH = 0.5,
        INIT = rep(0, d),
        VERBOSEFLAG = 0
    )
    tn

    cpp_ctrl <- cpp_ctrl_init
    cpp_ctrl$T_INIT <- cpp_ctrl$MAXT+1
    cpp_ctrl$MAXT <- 2*n
    cpp_ctrl$SAMPLING_WINDOW <- 1000
    cpp_ctrl$STEPSIZE <- tn$best_step
    cpp_ctrl$EACH = round(.25*n, 0)
    burn_theta <- tn$best_theta
    mod_obj <- fit_isingGraph3(
        DATA_LIST = list(DATA = as.matrix(mat), CONSTRAINTS = rep(T, d)),
        METHOD = 'recycle_hyper',
        CPP_CONTROL = cpp_ctrl,
        INIT = burn_theta, #par_init,
        VERBOSEFLAG = 0
    )

    path_theta <- get_tidy_path3(MOD_OBJ = mod_obj, PATH_LAB = 'path_av_theta') %>% filter(row_number()!=1)
    path_theta <- tn$best_path  %>% filter(iter != 1) %>% bind_rows(path_theta)
    path_hnll <- path_theta %>%
        mutate(
            hnll = map_dbl(path_av_theta, ~ising_graph_lik(.x))
        )
    tictoc::toc()
}

tol <- 0.001
path_hnll$diff <- c(0, path_hnll %>% pluck('hnll') %>% diff())
path_hnll %>%
    mutate(check = abs(diff)/hnll)
st_theta <- path_hnll %>%
    mutate(check = abs(diff)/hnll) %>%
    filter(iter != min(path_hnll$iter) & check < tol) %>%
    pluck('path_av_theta', 1)
stop_iter <- path_hnll %>%
    mutate(check = abs(diff)/hnll) %>%
    filter(iter != min(path_hnll$iter) & check < tol) %>%
    pluck('iter', 1)
path_hnll <- path_hnll %>% filter(iter <= stop_iter)
path_theta <- path_theta %>% filter(iter <= stop_iter)

gg_hnll <- path_hnll %>%
    filter(iter >=1) %>%
    ggplot(aes(x = iter/n, y = hnll)) +
    geom_rect(aes(xmin=-Inf,xmax=cpp_ctrl_init$BURN/n, ymin=-Inf,ymax=Inf),
              fill="lightgrey", color =NA) +
    geom_line() +
    #geom_vline(xintercept = opt_iter/n, linetype = 'dashed') +
    labs(x = 'Number of passes through the training partition', y = 'Holdout negative composite log-likelihood') +
    theme_bw()
gg_hnll

traj <- path_theta %>%
    mutate(path_av_theta = map(path_av_theta, ~ tibble(par = (1:length(.x)), val = .x))) %>%
    unnest('path_av_theta') %>%
    mutate(lab = if_else(par > p, 'Edge coefficients', 'Intercepts'), pass = iter/n)
#traj %>% pluck('pass') %>% unique()

gg_traj <- traj %>%
    ggplot(aes(x = pass, y = val, group = par)) +
    geom_rect(aes(xmin=-Inf,xmax=cpp_ctrl_init$BURN/n,ymin=-Inf,ymax=Inf),
              fill="lightgrey", color =NA) +
    geom_line(alpha = .8) +
    labs(x = 'Number of passes through the training partition', y = 'Stochastic estimates') +
    facet_grid(~lab) +
    theme_bw() +
    theme(text=element_text(size=15))
gg_traj
gg_checks <- ggpubr::ggarrange(gg_traj + labs(x = ' '),
                               gg_hnll + labs(x = ' '),
                               nrow = 1, widths = c(2,1) )
gg_checks <- ggpubr::annotate_figure(gg_checks, bottom = 'Number of passes through the training partition')
gg_checks
#ggsave(gg_checks, filename = 'gg_checks.pdf', width = 10, height = 4)

#### inference ####
system.time(
    invH <- sampleH(THETA = st_theta, DATA = as.matrix(mat), CONSTRAINTS = rep(T, d), INVERTFLAG = T)
)
diag(invH)
J <-  sampleJ(THETA = st_theta, DATA = as.matrix(mat), CONSTRAINTS = rep(T, d))
sandwich <- invH %*% J %*% invH

se <- sqrt(diag(invH/(stop_iter-mod_obj$control$BURN) + sandwich/n))
zscores <- st_theta/se

pval <- 2*pnorm(q = abs(zscores), lower.tail = F)
pval_holm <- p.adjust(p = pval, method = 'holm')

alpha <- .01;
theta_holm <- if_else(pval_holm < alpha, st_theta, 0)

#### Plots ####
library(ggraph)
library(tidygraph)
library(stringr)
library(jcolors)
library(ggiraph)

# Get the edges matrices
emat_holm <- stIsing::ising_from_theta_to_emat(theta_holm, p)

edges_tib <- expand_grid(
    from = 2:p,
    to = 1:(p-1)) %>%
    filter(from > to) %>%
    mutate(
        presence_holm = map2_dbl(from, to, ~emat_holm[.x, .y]!=0),
        intensity = map2_dbl(from, to, ~emat_holm[.x, .y]),
        sign = if_else(intensity > 0, 'positive', 'negative'))
edges_tib %>% pluck('presence_holm') %>% mean()
holm_tidy_net <- tbl_graph(
    nodes = dictionary,
    edges = edges_tib %>%
        filter(presence_holm!=0),
    directed = F )%>%
    mutate(centrality = centrality_degree(mode = 'in')) %>%
    filter(centrality>0)

colz <- jcolors::jcolors(palette = 'pal7')
names(colz) <- NULL
gg_chosen <- holm_tidy_net  %>%
    ggraph(layout = "kk") +
    geom_edge_link(aes(width=abs(intensity), linetype = sign), alpha = 0.2,     show.legend = F) +
    geom_node_point(aes(col = area), alpha = .5, size = 2) +
    geom_point_interactive(aes(x = x, y = y, size = (centrality), col = area, data_id = area, tooltip = txt), alpha =.8,  show.legend = FALSE) +
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_label(aes(label = label, col = area), repel = TRUE, alpha = .5, size = 3, label.size = .3, show.legend = F) +
    theme_graph(base_family = 'Helvetica', base_size = 15,   #text_colour = "#888B8D",
                plot_margin = margin(1, 1, 1, 1)
    )+
    scale_edge_linetype_manual(values = c('positive' = 'solid', 'negative' = 'dashed'))+
    labs(col= ' ')+
    theme(legend.position = 'right',
          #panel.background = element_rect('#888B8D', fill = NA),
          plot.subtitle=element_text(hjust=0.5, vjust = -.3)) +
    scale_color_manual(values = colz)

#ggsave(gg_chosen, filename = 'mental_health_struct.pdf', width = 10, height = 6)

grf <- girafe(ggobj = gg_chosen, width_svg = 10, height_svg = 10,
              options = list(opts_sizing(rescale = F)))
grf
#htmlwidgets::saveWidget(grf, "mental_health_struct.html", selfcontained = T)

# Estimation check
gg1 <- tibble(
    num = num_theta, # Numerical estimate computed on real data
    stoc = st_theta,  # Stochastic estimates computed on synthetic data generated from num_theta
    id = 1:d
    ) %>%
  filter(id > p) %>% 
    ggplot(aes(x = num, y = stoc)) +
    geom_abline(slope = 1, linetype = 'dashed')+
    geom_point(alpha = .8) +
  theme_bw() +
  labs( x = 'True parameters (Numerical estimates on original data)', y = 'Stochastic estimates on simulated data')

gg2 <- tibble(
  stoc = st_theta, # Numerical estimate computed on real data
  stoc_original = st_theta_original,  # Stochastic estimates computed on synthetic data generated from num_theta
  id = 1:d
) %>%
  filter(id > p) %>% 
  ggplot(aes(x = stoc_original, y = stoc)) +
  geom_abline(slope = 1, linetype = 'dashed')+
  geom_point(alpha = .8) +
  theme_bw() +
  labs( x = 'Stochastic estimates on original data', y = 'Stochastic estimates on simulated data')

gg_estimation_checks <- ggpubr::ggarrange(gg1, gg2)
#ggsave(gg_estimation_checks, filename = 'gg_estimation_checks.pdf', width = 10, height = 6)
