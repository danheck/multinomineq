#' Encompassing Bayes Factor
#'
#' Computes the Bayes factor for order-constrained product-binomial models
#' that are specified as a polytope via:  A*x <= b (see details).
#'
#' @inheritParams compute_cnml
#' @inheritParams select_nml
#' @param A a matrix with one row for each linear inequality constraint and one
#'   column for each of the binomial parameters. The parameter space is defined
#'   as all probabilities \code{x} that fulfill the order constraints  \code{A*x <= b}.
#' @param b a vector of the same length as the number of rows of \code{A}
#' @param M number of posterior samples drawn from the encompassing model
#' @param batch size of the batches into which computations are split to reduce memory load
#' @param steps integer vector that indicates at which rows the matrix \code{A} is split for a stepwise computation of the Bayes factor (see details). In this case, \code{M} can be a vector with the number of samples drawn in each step from the (partially) order-constrained models using Gibbs sampling
#' @param prior a vector with two positive numbers defining the shape parameters of the beta prior distributions for each binomial rate parameter
#' @details
#'
#'
#' The stepwise computation of the Bayes factor proceeds as follows: If the steps are defined as \code{steps=c(5,10)}, the BF is computed in three steps by comparing: (1) Unconstrained model vs. inequalities in \code{A[1:5,]}; (2)  inequalities in \code{A[1:5,]} vs \code{A[1:10,]}; (3)  inequalities in \code{A[1:10,]} vs. all inequalities in \code{A}.
#' @template ref_karabatsos2005
#' @template ref_regenwetter2014
#' @return a list with the Bayes factor (and log BF) of the order-constrained vs. encompassing (\code{bf_0e}) and the encompassing vs. order-constrained (\code{bf_e0}) model, and the number of posterior/prior samples that satisfied the order-constraints.
#' @examples
#' k <- c(0, 3, 2, 5, 3, 7)
#' n <- rep(10, 6)
#'
#' # linear order constraint:
#' # theta1 < theta2 < .... < .50
#' A <- matrix(c(1, -1, 0, 0, 0, 0,
#'               0, 1, -1, 0, 0, 0,
#'               0, 0, 1, -1, 0, 0,
#'               0, 0, 0, 1, -1, 0,
#'               0, 0, 0, 0, 1, -1,
#'               0, 0, 0, 0, 0, 1),
#'             ncol = 6, byrow = TRUE)
#' b <- c(0, 0, 0, 0, 0, .50)
#'
#' # check whether vectors are in polytope:
#' A %*% c(.05, .1, .12, .16, .19, .23) <= b
#' A %*% c(.05, .3, .12, .16, .19, .53) <= b
#'
#' # Bayes factor: unconstrained vs. constrained
#' compute_bf(k, n, A, b, c(1, 1), M = 2e5)
#' compute_bf(k, n, A, b, c(1, 1), M = 1e4, steps = c(2,4,5))
#' @export
compute_bf <- function(k, n, A, b, prior = c(1, 1),
                       M = 5e5, steps, batch = 5000){
  check_knAbprior(k, n, A, b, prior)

  if (missing(steps) || is.null(steps)){
    bfe <- as.list(encompassing_bf(k, n, A, b, prior, M, batch))
    res <- c(list("bf" = get_all_bfs(bfe$posterior, bfe$prior, M),
                  "se" = se_bf(bfe$posterior, bfe$prior, bfe$M)),
             bfe[2:4])
  } else {
    check_stepsA(steps, A)
    zeros <- rep(0, length(k))
    mpost <- encompassing_stepwise(k, n, A, b, prior, M, steps, batch)
    mprior <- encompassing_stepwise(zeros, zeros, A, b, prior, M, steps, batch)
    res <- list("bf" = get_all_bfs(mpost$counts, mprior$counts, M),
                "se" = se_bf(mpost$counts, mprior$counts, M),
                "posterior" = c(mpost$counts),
                "prior" = c(mprior$counts), "M" = M)
  }
  res
}

get_all_bfs <- function (post, prior, Mpost = 1, Mprior = Mpost){
  if (!is.null(dim(post)) && ncol(post) > 1){
    bf_0e <- apply(post / prior * matrix(Mprior / Mpost,
                                         nrow(post), ncol(post), byrow = TRUE),
                   1, prod)
  } else {
    bf_0e <- prod(post / prior * Mprior / Mpost)
  }
  c("bf_0e" = bf_0e, "bf_e0" = 1 / bf_0e,
    "log_bf_0e" = log(bf_0e), "log_bf_e0" = - log(bf_0e))
}
