#' Bayes Factor for Stochastic Dominance of Continuous Distributions
#'
#' Uses discrete bins (as in a histogram) to compute the Bayes factor in favor
#' of stochastic dominance of continuous distributions.
#'
#' @param x1 a vector with samples from the first random variable/experimental condition.
#' @param x2 a vector with samples from the second random variable/experimental condition.
#' @param breaks number of bins of histogram. See \code{\link[graphics]{hist}}.
#' @inheritParams stochdom_Ab
#' @param ... further arguments passed to \code{\link{bf_multinom}}. Note that the
#'    noninformative default prior \code{1/number_of_bins} is used.
#' @template ref_heathcote2010
#' @examples
#' x1 <- rnorm(300, 0, 1)
#' x2 <- rnorm(300, .5, 1) # dominates x1
#' x3 <- rnorm(300, 0, 1.2) # intersects x1
#'
#' plot(ecdf(x1))
#' lines(ecdf(x2), col = "red")
#' lines(ecdf(x3), col = "blue")
#'
#' b12 <- stochdom_bf(x1, x2, order = "<", M = 5e4)
#' b13 <- stochdom_bf(x1, x3, order = "<", M = 5e4)
#' b12$bf
#' b13$bf
#' @export
stochdom_bf <- function(x1, x2, breaks = "Sturges", order = "<", ...) {
  histo <- hist(c(x1, x2), breaks = breaks, plot = FALSE)
  b <- histo$breaks
  h1 <- hist(x1, breaks = b, plot = FALSE)
  h2 <- hist(x2, breaks = b, plot = FALSE)
  k1 <- h1$counts
  k2 <- h2$counts
  k <- c(k1, k2)
  B <- length(k1)
  bmat <- rbind("x1" = k1, "x2" = k2)
  colnames(bmat) <- paste0("b", 1:B)

  Ab <- stochdom_Ab(bins = B, conditions = 2, order = order)
  bf <- bf_multinom(k,
    options = rep(B, 2), A = Ab$A, b = Ab$b,
    prior = rep(1 / B, B * 2), ...
  )
  list(
    "breaks" = b, "k" = bmat,
    "A" = Ab$A, "b" = Ab$b,
    "bf" = bf
  )
}
