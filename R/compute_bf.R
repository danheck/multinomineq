#' Encompassing Bayes Factor
#'
#' Computes the Bayes factor for order-constrained product-binomial models
#' that are specified as a polytope via:  A*x <= b (see details).
#'
#' @inheritParams count_binomial
#' @details
#' For more control how many samples should be drawn from the prior/posterior, use \code{\link{count_binomial}}.
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

  pr <- count_binomial (A, b, 0, 0, prior, M, steps, batch)
  po <- count_binomial (A, b, k, n, prior, M, steps, batch)
  count_to_bf(po, pr)
}
