check_stoc_args <- function(ARGS, N, D){

    out <- ARGS

    if(is.null(ARGS$MAXT)) out$MAXT <- N
    if(is.null(ARGS$BURN)) out$BURN <- N*.25
    if(is.null(ARGS$STEPSIZE)) out$STEPSIZE <- 1
    if(is.null(ARGS$NU)) out$NU <- 1
    if(is.null(ARGS$SEED)) out$SEED <- 123
    if(is.null(ARGS$SCALEVEC)) out$SCALEVEC <- rep(1, D)
    if(is.null(ARGS$EACH)) out$EACH <- 1
    if(is.null(ARGS$HOLDOUTFLAG)) out$HOLDOUTFLAG <- F
    if(is.null(ARGS$PAR1)) out$PAR1 <- 1
    if(is.null(ARGS$PAR2)) out$PAR2 <- 1
    if(is.null(ARGS$PAR3)) out$PAR3 <- .501
    return(out)
}

#' Get parameters trajectories
#'
#' @param MOD_OBJ object from \link{fit_isingGraph3}
#' @param PATH_LAB Label of the path of interest
#'
#' @export
get_tidy_path3 <- function(MOD_OBJ, PATH_LAB){
    iters <- MOD_OBJ$fit$iter_idx
    path  <- Reduce(rbind, MOD_OBJ$fit[[PATH_LAB]])

    if(requireNamespace(c('dplyr'), quietly = TRUE)){
        if(PATH_LAB%in%c('path_nll')){
            out <- tibble::tibble(iter = iters) |>
                dplyr::mutate(
                    path_chosen = c(path,NA)
                )
        }else if(PATH_LAB%in%c('path_grad')){
            iters <- iters[-1]
            out <- tibble::tibble(iter = iters) |>
                dplyr::mutate(
                    path_chosen = split(t(path), rep(1:nrow(path), each = ncol(path)))
                )

        }else{
            out <- tibble::tibble(iter = iters) |>
                dplyr::mutate(
                    path_chosen = split(t(path), rep(1:nrow(path), each = ncol(path)))
                )
        }


        colnames(out) <- c('iter', PATH_LAB)

        return(out)
    }else{return(NULL)}

}

#' From vector to matrix
#'
#' @param par Parameter vector
#' @param p Number of nodes
#'
#' @export
ising_from_theta_to_emat <- function(par, p){
    emat <- matrix(0, p, p)
    counter <- p+1
    for (col in 1:(p-1)) {
        for (row in (col+1):p) {
            #cat('row:', row, 'col:', col, ' counter:', counter, '\n')
            emat[row, col] <- par[counter]
            counter <- counter + 1
        }
    }
    return(emat)
}

#' From matrix to vector
#' @param graph Graph matrix
#' @param intercepts Nodes intercepts
#'
#' @export
ising_from_graph_to_theta <- function(graph, intercepts){
    p <- length(intercepts)
    par <- intercepts

    counter <- p+1
    for (col in 1:(p-1)) {
        for (row in (col+1):p) {
            #cat('row:', row, 'col:', col, ' counter:', counter, '\n')
            par[counter] <- graph[row, col]
            counter <- counter + 1
        }
    }

    return(par)
}
