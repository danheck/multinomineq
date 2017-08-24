#' Number of Product-Binomial Prior/Posterior Samples in Polytope
#'
#' Draws prior/posterior samples for product-binomial data and counts how many samples are
#' inside the convex polytope defined by
#' (1) the inequalities A*x <= b or
#' (2) the convex hull over the vertices V.
#'
#' @inheritParams inside
#' @param k the number of Option B choices.
#'     The default \code{k=n=rep(0,ncol(A))} is equivalent to sampling from the prior.
#' @param n the number of choices per item type.
#' @param map optional: numeric vector of the same length as \code{k} with integers mapping the frequencies \code{k} to the free parameters/columns of \code{A}/\code{V}, thereby allowing for equality constraints (e.g., \code{map=c(1,1,2,2)}). Reversed probabilities \code{1-p} are coded by negative integers. Guessing probabilities of .50 are encoded by zeros. The default assumes different parameters for each item type: \code{map=1:ncol(A)}
#' @param M number of posterior samples drawn from the encompassing model
#' @param batch size of the batches into which computations are split to reduce memory load
#' @param steps integer vector that indicates at which rows the matrix \code{A} is split for a stepwise computation of the Bayes factor (see details). In this case, \code{M} can be a vector with the number of samples drawn in each step from the (partially) order-constrained models using Gibbs sampling
#' @param prior a vector with two positive numbers defining the shape parameters of the beta prior distributions for each binomial rate parameter.
#' @param start only if \code{steps} is defined: a vector with starting values
#'     within the polytope. If \code{steps = -1}, a random starting value is used
#'     (but a poitn within the polytope might not be found if the order constraints are strong).
#' @param progress whether a progress bar should be shown.
#'
#' @details
#' Useful to compute the encompassing Bayes factor for testing order constraints (see \code{\link{bf_binomial}}; Hojtink, 2011).
#'
#' The stepwise computation of the Bayes factor proceeds as follows:
#' If the steps are defined as \code{steps=c(5,10)}, the BF is computed in three steps by comparing:
#' (1) Unconstrained model vs. inequalities in \code{A[1:5,]};
#' (2) use posterior based on inequalities in \code{A[1:5,]} and check inequalities \code{A[6:10,]};
#' (3) sample from A[1:10,] and check inequalities in \code{A[11:nrow(A),]} (i.e., all inequalities).
#'
#' @return a list with the elements
#' \itemize{
#'     \item\code{integral}: estimated probability that samples are in polytope
#'     \item\code{count}: number of samples in polytope
#'     \item\code{M}: total number of samples
#'     \item\code{const_map}: logarithm of the binomial constants that have to be considered due to equality constraints
#' }
#' @examples
#' # linear order constraint:
#' # x1 < x2 < .... < x6 < .50
#' A <- matrix(c(1, -1, 0, 0, 0, 0,
#'               0, 1, -1, 0, 0, 0,
#'               0, 0, 1, -1, 0, 0,
#'               0, 0, 0, 1, -1, 0,
#'               0, 0, 0, 0, 1, -1,
#'               0, 0, 0, 0, 0, 1),
#'             ncol = 6, byrow = TRUE)
#' b <- c(0, 0, 0, 0, 0, .50)
#'
#' k <- c(0, 3, 2, 5, 3, 7)
#' n <- rep(10, 6)
#'
#' # check whether specific vector is in polytope:
#' A %*% c(.05, .1, .12, .16, .19, .23) <= b
#'
#' # count prior samples and compare to analytical result
#' prior <- count_binomial(0, 0, A, b, M = 5e5)
#' prior
#' (.50)^6 / factorial(6)
#'
#' # count posterior samples and get Bayes factor
#' posterior <- count_binomial(k, n, A, b,
#'                             M=c(1e5, 2e4), steps = 2)
#' posterior$integral / prior$integral  # Bayes factor
#' count_to_bf(posterior, prior)
#' @template ref_hoijtink2011
#' @template ref_fukuda2004
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
count_binomial <- function (k, n, A, b, V, map, prior = c(1, 1),
                            M = 10000, steps, batch = 10000,
                            start = -1, progress = TRUE){

  if (missing(A)) A <- V
  aggr <- map_k_to_A(k, n, A, map, prior)
  k <- aggr$k
  n <- aggr$n

  check_Mbatch(M, batch)

  if (!missing(b)){
    check_Abknprior(A, b, k, n, prior)
    if (missing(steps) || is.null(steps) || length(steps) == 0){
      cnt <- as.list(count_binomial_cpp(k, n, A, b, prior, M, batch, progress))
    } else {
      check_stepsA(steps, A)
      cnt <- count_stepwise(k, n, A, b, prior, M, steps, batch, start, progress)
    }

  } else if (!missing(V)){
    count <- 0
    m <- M
    a <- c(rbind(k + prior[1], n - k + prior[2]))
    while (m > 0 ){
      X <- rpdirichlet_free(m, a, rep(2, ncol(V)))
      count <- count + sum(inside_V(X, V))
      m <- m - batch
    }
    cnt <- list("integral" = count/M, "count" = count, "M" = M,
                "const_map_0e" = aggr$const_map_0e)
  } else {
    stop("A/b or V must be provided.")
  }
  cnt$const_map_0e <- aggr$const_map_0e
  cnt
}

# count_V <- function(X, V){
#   apply(X, 1, inside_V, V = V)
# }


