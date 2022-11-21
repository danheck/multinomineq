#' Random Generation for Independent Multinomial Distributions
#'
#' Generates random draws from independent multinomial distributions (= product-multinomial \code{pmultinom}).
#'
#' @param prob vector with probabilities or a matrix with one probability vector per row.
#'    For \code{rpbinom}: probabilities of a success for each option.
#'    For \code{rpmultinom}: probabilities of all categories excluding
#'    the last category for each option (cf. \code{drop_fixed}).
#'    See also \code{\link{sampling_binom}} and \code{\link{sampling_multinom}}.
#' @param n integer vector, specifying the number of trials for each binomial/multinomial distribution
#'    Note that this is the \code{size} argument in \code{rmultinom}, cf. \code{\link[stats]{Multinom}}.
#' @param drop_fixed whether the output matrix includes the last probability for each category
#'     (which is not a free parameter since probabilities must sum to one).
#' @inheritParams count_multinom

#' @return a matrix with one vector of frequencies per row. For \code{rpbinom}, only
#'    the frequencies of 'successes' are returned, whereas for \code{rpmultinom}, the
#'    complete frequency vectors (which sum to \code{n} within each option) are returned.
#' @examples
#' # 3 binomials
#' rpbinom(prob = c(.2, .7, .9), n = c(10, 50, 30))
#'
#' # 2 and 3 options:  [a1,a2,  b1,b2,b3]
#' rpmultinom(
#'   prob = c(a1 = .5, b1 = .3, b2 = .6),
#'   n = c(10, 20), options = c(2, 3)
#' )
#' # or:
#' rpmultinom(
#'   prob = c(a1 = .5, a2 = .5, b1 = .3, b2 = .6, b3 = .1),
#'   n = c(10, 20), options = c(2, 3),
#'   drop_fixed = FALSE
#' )
#'
#' # matrix with one probability vector per row:
#' p <- rpdirichlet(
#'   n = 6, alpha = c(1, 1, 1, 1, 1),
#'   options = c(2, 3)
#' )
#' rpmultinom(prob = p, n = c(20, 50), options = c(2, 3))
#' @export
rpbinom <- function(prob, n) {
  k <- rpmultinom(
    prob = prob, n = n,
    options = rep(2, length(prob)), drop_fixed = TRUE
  )
  drop_fixed(k, options = 2)
}

#' @rdname rpbinom
#' @export
rpmultinom <- function(prob, n, options, drop_fixed = TRUE) {
  if (length(n) == 1) {
    n <- rep(n, sum(options))
  }
  if (is.null(dim(prob))) {
    prob <- t(prob)
  }
  check_probko(
    prob = prob, k = rep(0, sum(options)),
    options = options, drop_fixed = drop_fixed
  )

  if (drop_fixed) prob <- add_fixed(prob, options = options)
  rr <- rpm_mat(prob, n, options)
  colnames(rr) <- colnames(prob)
  if (is.null(colnames(rr))) {
    colnames(rr) <- index_mult(options, fixed = TRUE)
  }
  rr
}
