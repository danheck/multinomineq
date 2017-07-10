#' Posterior Samples for Ordered-Constrained Models of Binary Choice
#'
#' Uses Gibbs sampling to draw posterior samples for order-constrained models of binary choice.
#'
#' @inheritParams compute_marginal
#' @inheritParams compute_bf
#' @param M number of posterior samples
#' @param start starting vectors. Must be in the convex polytope (that is, \code{A*start <= b})
#'
#' @details
#' Can be used to obtain posterior samples for binary random utility models that assume a mixture over predefined
#' preference orders/vertices that jointly define a convex polytope via (A * x <= b).
#' @return a matrix with posterior samples for the probabilities to choose Option B for each item type
#' @template ref_myung2005
#' @examples
#' A <- matrix(c(1, 0, 0,   # x1 < .50
#'               1, 1, 1,   # x1+x2+x3 < 1
#'               0, 2, 2,   # 2*x2+2*x3 < 1
#'               0, -1, 0,  # x2 > .2
#'               0, 0, 1),  # x3 < .1
#'             ncol = 3, byrow = TRUE)
#' b <- c(.5, 1, 1, -.2, .1)
#' samp <- sampling_binary(c(5,12,7), c(20,20,20), A, b)
#' head(samp)
#' colMeans(samp)
#' apply(samp, 2, plot, type = "l", ylim = c(0, 1))
#' @export
sampling_binary <- function (k, n, A, b, prior = c(1, 1), M = 5000, start){
  check_knAbprior(k, n, A, b, prior)
  if (missing(start)){
    start <- -1  # TODO: find start value
  } else if (length(start) != ncol(A) || !all(A %*% start <= b)) {
    stop ("'start' must be in the convex polytope:  A*start <= b")
  }
  samples <- sampling_binary_cpp(k, n, A, b, prior, M, start)
  samples
}
