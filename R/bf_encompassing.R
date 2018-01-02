#' Bayes Factor for Linear Order Constraints
#'
#' Computes the Bayes factor for linear order-constraints (specified via: A*x <= b)
#' for product-binomial/-multinomial models.
#'
#' @param ... further arguments passed to \code{\link{count_binom}}/\code{\link{count_multinom}}
#'    (e.g., \code{M}, \code{steps}).
#' @inheritParams count_binom
#' @inheritParams count_to_bf
#' @details
#' For more control, use \code{\link{count_binom}} to specifiy how many samples
#' should be drawn from the prior and posterior, respectively.
#' This is especially recommended if the same prior distribution (and thus the same prior probability/integral)
#' is used for computing BFs for multiple data sets that differ only in \code{k} and \code{n}.
#' In this case, the prior probability/proportion of the parameter space in line with the
#' inequality constraints can be computed once with high precision (or even analytically),
#' and only the posterior probability/proportion needs to be estimated separately for each unique vector \code{k}.
#'
#' @template ref_karabatsos2005
#' @template ref_regenwetter2014
#' @template return_bf
#' @seealso \code{\link{count_binom}} and \code{\link{count_multinom}} for
#'     for more control on the number of prior/posterior samples and
#'     \code{\link{bf_nonlinear}} for nonlinear order constraints.
#' @examples
#' k <- c(0, 3, 2, 5, 3, 7)
#' n <- rep(10, 6)
#'
#' # linear order constraint:
#' # prob_1 < prob_2 < .... < .50
#' A <- matrix(c(1, -1, 0, 0, 0, 0,
#'               0, 1, -1, 0, 0, 0,
#'               0, 0, 1, -1, 0, 0,
#'               0, 0, 0, 1, -1, 0,
#'               0, 0, 0, 0, 1, -1,
#'               0, 0, 0, 0, 0, 1),
#'             ncol = 6, byrow = TRUE)
#' b <- c(0, 0, 0, 0, 0, .50)
#'
#' # Bayes factor: unconstrained vs. constrained
#' bf_binom(k, n, A, b, prior=c(1, 1), M=5e4)
#' bf_binom(k, n, A, b, prior=c(1, 1), M=5000, steps=c(2,4,5))
#' bf_binom(k, n, A, b, prior=c(1, 1), M=1000, cmin = 500)
#' @export
bf_binom <- function(k, n, A, b, V, map, prior = c(1, 1), log = FALSE, ...){
  pr <- count_binom(0, 0, A, b, V, map = map, prior = prior, ...)
  po <- count_binom(k, n, A, b, V, map = map, prior = prior, ...)
  count_to_bf(po, pr, log = log)
}

#' @inheritParams count_multinom
#' @rdname bf_binom
#' @export
bf_multinom <- function(k, options, A, b, V, prior = rep(1, sum(options)),
                        log = FALSE, ...){
  pr <- count_multinom(0, options, A, b, V, prior, ...)
  po <- count_multinom(k, options, A, b, V, prior, ...)
  count_to_bf(po, pr, log = log)
}
