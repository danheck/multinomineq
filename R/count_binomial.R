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
#' @param steps an integer vector that indicates at which rows the matrix \code{A}
#'     is split for a stepwise computation of the Bayes factor (see details).
#'     In this case, \code{M} can be chosen to be a vector with the number of samples drawn
#'     in each step from the (partially) order-constrained models  using Gibbs sampling.
#'     If \code{cmin>0}, steps are chosen by default as: \code{steps=1:nrow(A)}.
#' @param prior a vector with two positive numbers defining the shape parameters of the beta prior distributions for each binomial rate parameter.
#' @param start only if \code{steps} is defined: a vector with starting values
#'     within the polytope. If \code{steps = -1}, a random starting value is used
#'     (but a point within the polytope might not be found if the order constraints are strong).
#' @param cmin if \code{cmin>0}: minimum number of counts per step in the automatic stepwise procedure.
#' @param maxiter if \code{steps="auto"}: maximum number of sampling runs in the automatic stepwise procedure.
#' @param progress whether a progress bar should be shown.
#'
#' @details
#' Useful to compute the encompassing Bayes factor for testing order constraints (see \code{\link{bf_binom}}; Hojtink, 2011).
#'
#' The stepwise computation of the Bayes factor proceeds as follows:
#' If the steps are defined as \code{steps=c(5,10)}, the BF is computed in three steps by comparing:
#' (1) Unconstrained model vs. inequalities in \code{A[1:5,]};
#' (2) use posterior based on inequalities in \code{A[1:5,]} and check inequalities \code{A[6:10,]};
#' (3) sample from A[1:10,] and check inequalities in \code{A[11:nrow(A),]} (i.e., all inequalities).
#'
#' @return a matrix with the columns
#' \itemize{
#'     \item\code{count}: number of samples in polytope / that satisfy order constraints
#'     \item\code{M}: the  total number of samples in each step
#'     \item\code{steps}: the \code{"steps"} used to sample from the polytope
#'         (i.e., the row numbers of \code{A} that were considered  stepwise)
#' }
#' with the attributes:
#' \itemize{
#'    \item\code{integral}: estimated probability that samples are in polytope
#'    \item\code{const_map}: logarithm of the binomial constants that
#'           have to be considered due to equality constraints
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
#' prior <- count_binom(0, 0, A, b, M = 1e4, steps = 1:4)
#' prior  # volume = attribute "integral"
#' (.50)^6 / factorial(6)
#'
#' # count posterior samples stepwise
#' posterior <- count_binom(k, n, A, b, M=1e4, steps=1:4)
#' count_to_bf(posterior, prior)
#'
#' # automatic stepwise algorithm
#' prior <- count_binom(0, 0, A, b, cmin = 1000)
#' posterior <- count_binom(k, n, A, b, cmin = 1000)
#' count_to_bf(posterior, prior)

#' @template ref_hoijtink2011
#' @template ref_fukuda2004
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
count_binom <- function (k, n, A, b, V, map, prior = c(1, 1),
                         M = 10000, steps, start = -1,
                         cmin = 0, maxiter = 500, progress = TRUE){

  check_Mminmax(M, cmin, maxiter)
  if (missing(A)) A <- V
  if (start[1] == -1) start <- find_inside(A, b)

  aggr <- map_k_to_A(k, n, A, map, prior)
  k <- aggr$k
  n <- aggr$n


  if (!missing(b)){
    check_Abknprior(A, b, k, n, prior)
    if (cmin == 0 && missing(steps)){
      count <- count_bin(k, n, A, b, prior, M, batch = BATCH, progress)
    } else if (cmin > 0){
      if (missing(steps)) steps <- seq(1, nrow(A))
      else steps <- check_stepsA(steps, A)
      zeros <- rep(0, length(steps))
      count <- count_auto_bin(k, n, A, b, prior, zeros, zeros, steps,
                              M_iter = M, cmin = cmin, maxiter = maxiter, start, progress)
    } else {
      steps <- check_stepsA(steps, A)
      count <- count_stepwise_bin(k, n, A, b, prior, M, steps, batch = BATCH, start, progress)
    }

  } else if (!missing(V)){
    count <- 0
    m <- M
    a <- c(rbind(k + prior[1], n - k + prior[2]))
    while (m > 0 ){
      X <- rpdirichlet_free(m, a, rep(2, ncol(V)))
      count <- count + sum(inside_V(X, V))
      m <- m - BATCH
    }
    count <- cbind("count" = count, "M" = M, "steps" = NA)
  } else {
    stop("A/b or V must be provided.")
  }
  attr(count, "integral") <- get_integral(count)
  attr(count, "const_map_0e") <- aggr$const_map_0e
  count
}

get_integral <- function(count){
  prod(count[,"count"]/count[,"M"])
}
