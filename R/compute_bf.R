#' Encompassing Bayes Factor
#'
#' Computes the Bayes factor for order-constrained product-binomial models
#' that are specified as a polytope via:  A*x < b (see details).
#'
#' @inheritParams compute_cnml
#' @inheritParams select_nml
#' @param A a matrix with one row for each linear inequality constraint and one
#'   column for each of the binomial parameters.
#' @param b a vector of the same length as the number of rows of \code{A}
#' @param M number of posterior samples drawn from the encompassing model
#' @param batch size of the batches into which computations are split to reduce memory load
#' @param steps integer vector that indicates at which rows the matrix \code{A} is split for a stepwise computation of the Bayes factor (see details).
#'
#' @details
#'
#'
#' The stepwise computation of the Bayes factor proceeds as follows: If the steps are defined as \code{steps=c(5,10)}, the BF is computed in three steps by comparing: (1) Unconstrained model vs. inequalities in \code{A[1:5,]}; (2)  inequalities in \code{A[1:5,]} vs \code{A[1:10,]}; (3)  inequalities in \code{A[1:10,]} vs. all inequalities in \code{A}.
#' @template ref_karabatsos2005
#' @template ref_regenwetter2014
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
#' compute_bf(k, n, A, b, c(1, 1))
#' @export
compute_bf <- function(k, n, A, b, prior = c(1, 1),
                       M = 2e5, steps, batch = 5000){
  check_knpcp(k, n, k, .5, prior)
  check_kAb(k, A, b)

  if (missing(steps) || is.null(steps)){
    res <- as.list(encompassing_bf(k, n, A, b, prior, M, batch))
    se <- se_bf(res$posterior, res$prior, res$M)
    res <- c(res[1], list("se" = se), res[2:4])
  } else {
    check_stepsA(steps, A)
    zeros <- rep(0, length(k))
    m1 <- encompassing_stepwise(k, n, A, b, prior, M, steps, batch)
    m0 <- encompassing_stepwise(zeros, zeros, A, b, prior, M, steps, batch)
    se <- se_bf(m1$counts, m0$counts, M)
    res <- list("bf" = prod(m1$counts/m0$counts), "se" = se,
                "posterior" = c(m1$counts), "prior" = c(m0$counts), "M" = M)
  }
  res
}
