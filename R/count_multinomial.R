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
#' @inheritParams count_binom
#' @return a list with the elements
#' @return a matrix with the columns
#' \itemize{
#'     \item\code{count}: number of samples in polytope / that satisfy order constraints
#'     \item\code{M}: the  total number of samples in each step
#'     \item\code{steps}: the \code{"steps"} used to sample from the polytope
#'         (i.e., the row numbers of \code{A} that were considered stepwise)
#' }
#' with the attributes:
#' \itemize{
#'    \item\code{integral}: estimated probability that samples are in polytope
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
#' prior <- count_multinom(0, options, A, b, M = 2e4)
#' posterior <- count_multinom(k, options, A, b, M = 2e4)
#' count_to_bf(posterior, prior)
#'
#' bf_multinom(k, options, A, b, M=1e5)
#' bf_multinom(k, options, A, b, steps = 2, cmin = 5000)
#' @export
count_multinom <- function (k = 0, options, A, b, V, prior = rep(1, sum(options)),
                            M = 5000, steps, start = -1,
                            cmin = 0, maxiter = 500, progress = TRUE){
  if (length(k) == 1 && k == 0) k <- rep(0, sum(options))
  if (start[1] == -1) start <- find_inside(A, b)

  check_Abokprior(A, b, options, k, prior)
  check_Mminmax(M, cmin, maxiter)

  if (cmin == 0 && missing(steps)){
    count <- count_mult(k, options, A, b, prior, M, batch = BATCH, progress)
  } else {
    tmp <- Ab_multinom(options, A, b)  # sum-to-1 constraints (for Gibbs sampler!)
    A <- tmp$A
    b <- tmp$b
    if (missing(steps)) steps <- seq(1, nrow(A))
    steps <- check_stepsA(steps, A)

    if (cmin > 0){
      zeros <- rep(0, length(steps))
      count <- count_auto_mult(k, options, A, b, prior, zeros, zeros, steps,
                               M_iter = M, cmin = cmin, maxiter = maxiter, start, progress)
    } else {
      count <- count_stepwise_multi(k, options, A, b, prior, M,
                                    steps, batch = BATCH, start, progress)
    }
  }
  attr(count, "integral") <- get_integral(count)
  count
}
