#' Random Sample for Product-Multinomial Distribution
#'
#' Generate random samples from independent multinomial distributions.
#'
# ' @param M number of random vectors to draw.
#' @param prob vector with probability parameters or a matrix with one vector per row,
#'    (in which case \code{M} is ignored).
#'    Values must be nonnegative and smaller than 1 for \code{rpbinom} and
#'    nonnegative and sum to one within each option for \code{rpmultinom}.
#'    See \code{\link{sampling_binom}} and \code{\link{sampling_multinom}}.
#' @param n integer vector, specifying the number of trials for each binomial/multinomial distribution
#'    (note: this is the \code{size} argument in \code{\link[stats]{rmultinom}}).
#' @param options number of parameters/categories of each multinomial distribution.
#' @return a matrix with one vector of frequencies per row. For \code{rpbinom}, only
#'    the frequencies of 'successes' are returned, whereas for \code{rpmultinom}, the
#'    complete frequency vectors (which sum to \code{n} within each option) are returned.
#' @examples
#' # 5 binomials
#' rpbinom(c(.2, .7, .9), c(10, 50, 30))
#'
#' # 2 and 3 options:  [a1,a2,  b1,b2,b3]
#' options <- c(2, 3)
#' rpmultinom(c(.5,.5,   .3,.6,.1), c(100, 200), options)
#'
#' prob <- rpdirichlet(5, rep(1,5), options)
#' rpmultinom(prob, c(20, 50), options)
#' @export
rpbinom <- function(prob, n){
  if (length(n) == 1) n <- rep(n, length(prob))
  if (is.null(dim(prob)))
    prob <- matrix(prob, 1, length(prob), byrow = TRUE)
  if (is.null(dim(n)))
    n <- matrix(n, 1, length(n), byrow = TRUE)

  matrix(rbinom(length(prob), n, prob), nrow = nrow(prob),
         dimnames = list(NULL, colnames(prob)))
}

#' @rdname rpbinom
#' @export
rpmultinom <- function(prob, n, options){
  if (length(n) == 1) n <- rep(n, sum(options))
  if (is.null(dim(prob)))
    prob <- matrix(prob, 1, length(prob), byrow = TRUE)

  check_prob(prob)
  check_o(options)
  rr <- rpm_mat(prob, n, options)
  colnames(rr) <- colnames(prob)
  rr
}
