#' Number of Product-Multinomial Prior/Posterior Samples in Polytope
#'
#' Draws prior/posterior samples for product-multinomial data and counts how many samples are
#' inside the convex polytope defined by
#' (1) the inequalities A*x <= b or
#' (2) the convex hull over the vertices V.
#'
#' @param A a matrix defining the convex polytope via A*x <= b.
#'    The columns of A do not include the last choice option per item type and
#'    thus the number of columns must be equal to \code{sum(options-1)} (e.g., the column order of \code{A}
#'    for \code{k = c(a1,a2,a2, b1,b2)} is \code{c(a1,a2, b1)}).
#' @param options the number of choice options per item type, e.g., \code{c(3,2)} for a ternary and binary item.
#'     The sum of \code{options} must be equal to the length of \code{k}.
#' @param k the number of choices for each alternative ordered by item type, e.g.
#'     \code{c(a1,a2,a3, b1,b2)}.
#'     The default \code{k=0} is equivalent to sampling from the prior.
#' @param prior the prior parameters of the Dirichlet-shape parameters.
#'    Must have the same length as \code{k}.
#' @inheritParams inside
#' @return a list with the elements
#' \itemize{
#'     \item\code{integral}: estimated probability that samples are in polytope
#'     \item\code{count}: number of samples in polytope
#'     \item\code{M}: total number of samples
#' }
#' @template ref_hoijtink2011
#' @examples
#' ### frequencies:
#' #           (a1,a2,a3, b1,b2,b3,b4)
#' k <-       c(1,5,9,    5,3,7,8)
#' options <- c(3,        4)
#'
#' ### linear order constraints
#' # a1<a2<a3   AND   b2<b3<.50
#' # (note: a2<a3 <=> a2<1-a1-a2 <=> a1+2*a2<1)
#' # matrix A:
#' #             (a1,a2, b1,b2,b3)
#' A <- matrix(c(1, -1, 0,  0,  0,
#'               1,  2, 0,  0,  0,
#'               0,  0, 0,  1, -1,
#'               0,  0, 0,  0,  1),
#'             ncol = sum(options-1), byrow = TRUE)
#' b <- c(0, 1, 0, .50)
#'
#' # count prior and posterior samples and get BF
#' prior <- count_multinomial(0, options, A, b, M = 2e5)
#' posterior <- count_multinomial(k, options, A, b, M = c(2e5))
#' posterior$integral / prior$integral  # BF for constraints
#' count_to_bf(posterior, prior)
#'
#' bf_multinomial(k, options, A, b, M=1000)
#' @export
count_multinomial <- function (k = 0, options, A, b, V,
                               prior = rep(1, sum(options)),
                               M, steps, batch = 10000){
  if (length(k) == 1 && k == 0)
    k <- rep(0, sum(options))
  check_Abokprior(A, b, options, k, prior)
  check_Mbatch(M, batch)

  if (missing(steps) || is.null(steps) || length(steps) == 0){
    cnt <- as.list(count_multinomial_cpp(k, options, A, b, prior, M, batch))
  } else {
    check_stepsA(steps, A)
    cnt <- count_stepwise_multi(k, n, A, b, prior, M, steps, batch)
  }
  cnt
}
