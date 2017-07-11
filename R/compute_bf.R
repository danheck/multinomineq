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
#'
#' @details
#' The stepwise computation of the Bayes factor proceeds as follows: If the steps are defined as \code{steps=c(5,10)}, the BF is computed in three steps by comparing: (1) Unconstrained model vs. inequalities in \code{A[1:5,]}; (2) use posterior based on inequalities in \code{A[1:5,]} and check inequalities \code{A[6:10,]}; (3) sample from A[1:10,] and check inequalities in \code{A[11:nrow(A),]} (i.e., all inequalities).
#'
#' For more control how many samples should be drawn from the prior/posterior, use \code{\link{count_polytope}}.
#' This is especially recommended if the same prior distribution (and thus the same integral)
#' is used for computing BFs for multiple data sets that differ only in \code{k} and \code{n}.
#'
#' @template ref_karabatsos2005
#' @template ref_regenwetter2014
#' @template return_bf
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
                       M = 5e5, steps, batch = 10000){

  pr <- count_polytope (A, b, 0, 0, prior, M, steps, batch)
  po <- count_polytope (A, b, k, n, prior, M, steps, batch)
  count_to_bf(po, pr)
}
