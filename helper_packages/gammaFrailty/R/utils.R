rmvn <- function(SAMPLE_SIZE, VAR){
    dim <- ncol(VAR)
    sample <- t(t(chol(VAR))%*%matrix(rnorm(dim*SAMPLE_SIZE), dim, SAMPLE_SIZE))

    return(sample)
}

#' theta reparameterisation
#'
#' From constrained  to unconstrained parameter vector
#'
#' @param par Constrained parameter vector
#'
#' @export
partorepar <- function(par){
    repar <- par
    repar[1] <- -log( repar[1])
    #repar[1] <- log( repar[1])
    #repar[1] <- zofr( repar[1]-1.1)
    repar[2] <- zofr_cpp(repar[2])
    return(repar)
}

#' theta reparameterisation
#'
#' From unconstrained  to constrained parameter vector
#'
#' @param repar Unconstrained parameter vector
#'
#' @export
repartopar <- function(repar){
    par <- repar
    par[1] <- exp(-par[1])
    par[2] <- rofz_cpp(par[2])
    return(par)
}

check_SCSD_args <- function(ARGS, N, D){

    out <- ARGS

    if(is.null(ARGS$MAXT)) out$MAXT <- round(N,0)
    if(is.null(ARGS$BURN)) out$BURN <- 0
    if(is.null(ARGS$STEPSIZE)) out$STEPSIZE <- 1
    if(is.null(ARGS$NU)) out$NU <- 1
    if(is.null(ARGS$SEED)) out$SEED <- 123
    if(is.null(ARGS$SCALEVEC)) out$SCALEVEC <- rep(1, D)
    if(is.null(ARGS$PAR1)) out$PAR1 <- 1
    if(is.null(ARGS$PAR2)) out$PAR2 <- 1
    if(is.null(ARGS$PAR3)) out$PAR3 <- .501


    return(out)
}


#' Get parameters trajectories
#'
#' @param MOD_OBJ object from \link{fit_gammaFrailty2}
#' @param PATH_LAB Label of the path of interest
#' @param USEREPARTOPAR TRUE to reparametrize estimate to constrained parameterization.
#'
#' @export
get_tidy_path <- function(MOD_OBJ, PATH_LAB, USEREPARTOPAR){
    if(requireNamespace(c('dplyr','purrr'), quietly = TRUE)){
        iters <- MOD_OBJ$iterations_subset
        path  <- MOD_OBJ$fit[[PATH_LAB]]

        if(PATH_LAB%in%c('path_nll')){
            out <- dplyr::tibble(iter = iters)  |>
                dplyr::mutate(
                    path_chosen = c(path,NA)
                )
        }else if(PATH_LAB%in%c('path_grad')){
            out <- dplyr::tibble(iter = iters)  |>
                dplyr::mutate(
                    path_chosen = split(t(rbind(path, rep(NA, ncol(path)))), rep(1:(nrow(path)+1), each = ncol(path)))
                )

        }else{
            out <- dplyr::tibble(iter = iters)  |>
                dplyr::mutate(
                    path_chosen = split(t(path), rep(1:nrow(path), each = ncol(path)))
                )
        }


        if(USEREPARTOPAR){
            out <- out  |>
                dplyr::mutate(path_chosen = purrr::map(path_chosen, ~repartopar(.x)))
        }

        colnames(out) <- c('iter', PATH_LAB)

        return(out)
    }else{return(NULL)}

}

# sample_weights <- function(PROB, N, TRIPLETS, SEED, ITER){
#     set.seed(SEED + ITER)
#     sel <- NULL
#     while (purrr::is_empty(sel)) {
#         sel <- which(rbinom(nrow(TRIPLETS), 1, prob = PROB)==1)
#         out <- TRIPLETS[sel,]
#         if(is.vector(out)) out <- matrix(out, 1, 3)
#         out <- cbind(out, t = ITER)
#     }
#
#     return(out)
# }

# sample_optimisation_weights <- function(N, P, MAXITER, PAIRS_RANGE = 100, SEED = 123){
#     k_ref <- expand_grid(j = 0:(P-1), jp = 0:(P-1))  |>
#         filter(((jp - j) <= PAIRS_RANGE) & ((jp - j)>0))
#     K <- nrow(k_ref)
#     message(paste0('Dimension of pairs pool: ', K, ' x ', N))
#     triplets_ref <- k_ref  |>
#         expand_grid(i=0:(N-1))  |>
#         select(i, j, jp)  |>
#         as.matrix()
#
#     w_list <- purrr::map(1:MAXITER, ~sample_weight(PROB = 1/n, N = n, TRIPLETS = triplets_ref, SEED = SEED, ITER = .x))
#     out <- purrr::reduce(w_list, rbind)
#
#     return(out)
# }
