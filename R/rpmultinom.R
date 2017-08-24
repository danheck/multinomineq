#' Random Sample for Product-Multinomial Distribution
#'
#' Generate random samples from independent multinomial distributions.
#' @param theta vector with probability parameters. Must sum to one within each option.
#'    Can be a matrix with one parameter vector per row.
#'    See \code{\link{sampling_binomial}} and \code{\link{sampling_multinomial}}.
#' @param n integer vector, specifying the total number of objects for each multinomial distribution.
#' @param options number of parameters/categories of each multinomial distribution.
#' @examples
#' # 5 binomials
#' rpbinom(c(.2, .7, .9, .1, .5), rep(20, 5))
#'
#' # 2 and 3 options:  [a1,a2,  b1,b2,b3]
#' rpmultinom(c(.5,.5, .3,.6,.1), c(100, 200), c(2, 3))
#'
#' theta <- rpdirichlet(5, rep(1,5), c(2,3))
#' rpmultinom(theta, c(20, 20), c(2, 3))
#' @export
rpmultinom <- function(theta, n, options){
  if(is.null(dim(theta)))
    theta <- matrix(theta, 1)
  check_theta(theta)
  check_o(options)
  rpm_mat(theta, n, options)
}

#' @rdname rpmultinom
#' @export
rpbinom <- function(theta, n){
  options <- rep(2, length(theta))
  if (length(n) == 1)
    n <- rep(n, length(theta))
  if(is.null(dim(theta)))
    theta <- matrix(theta, 1)
  k <- rpmultinom(free_to_full(theta, options), n, options)
  drop_fixed(k)
}
