#' Count How Many Samples Satisfy Linear Inequalities (Binomial)
#'
#' Draws prior/posterior samples for product-binomial data and counts how many samples are
#' inside the convex polytope defined by
#' (1) the inequalities \code{A*x <= b} or
#' (2) the convex hull over the vertices \code{V}.
#'
#' @inheritParams inside
#' @param k vector of observed response frequencies.
#' @param n the number of choices per item type.
#'     If \code{k=n=0}, Bayesian inference is relies on the prior distribution only.
#' @param map optional: numeric vector of the same length as \code{k} with integers
#'     mapping the frequencies \code{k} to the free parameters/columns of \code{A}/\code{V},
#'     thereby allowing for equality constraints (e.g., \code{map=c(1,1,2,2)}).
#'     Reversed probabilities \code{1-p} are coded by negative integers.
#'     Guessing probabilities of .50 are encoded by zeros. The default assumes
#'     different parameters for each item type: \code{map=1:ncol(A)}
#' @param M number of posterior samples drawn from the encompassing model
#' @param steps an integer vector that indicates the row numbers at which the matrix \code{A}
#'     is split for a stepwise computation of the Bayes factor (see details).
#'     \code{M} can be a vector with the number of samples drawn
#'     in each step from the (partially) order-constrained models  using Gibbs sampling.
#'     If \code{cmin>0}, samples are drawn for each step until \code{count[i]>=cmin}.
#' @param prior a vector with two positive numbers defining the shape parameters
#'     of the beta prior distributions for each binomial rate parameter.
#' @param start only relevant if \code{steps} is defined or \code{cmin>0}:
#'     a vector with starting values in the interior of the polytope.
#'     If missing, an approximate maximum-likelihood estimate is used.
#' @param cmin if \code{cmin>0}: minimum number of counts per step in the automatic stepwise procedure.
#'     If \code{steps} is not defined, \code{steps=c(1,2,3,4,...)} by default.
#' @param maxiter if \code{cmin>0}: maximum number of sampling runs in the automatic stepwise procedure.
#' @param burnin number of burnin samples per step that are discarded. Since the
#'     maximum-likelihood estimate is used as a start value (which is updated for each step in
#'     the stepwise procedure in \code{count_multinom}), the \code{burnin}
#'     number can be smaller than in other MCMC applications.
#' @param progress whether a progress bar should be shown (if \code{cpu=1}).
#' @param cpu either the number of CPUs used for parallel sampling, or a parallel
#'     cluster  (e.g., \code{cl <- parallel::makeCluster(3)}).
#'     All arguments of the function call are passed directly to each core,
#'     and thus the total number of samples is \code{M*number_cpu}.
#'
#' @details
#' Counts the number of random samples drawn from beta distributions that satisfy
#' the convex, linear-inequalitiy constraints. The function is useful to compute
#' the encompassing Bayes factor for testing inequality-constrained models
#' (see \code{\link{bf_binom}}; Hojtink, 2011).
#'
#' The stepwise computation of the Bayes factor proceeds as follows:
#' If the steps are defined as \code{steps=c(5,10)}, the BF is computed in three steps by comparing:
#' (1) Unconstrained model vs. inequalities in \code{A[1:5,]};
#' (2) use posterior based on inequalities in \code{A[1:5,]} and check inequalities \code{A[6:10,]};
#' (3) sample from A[1:10,] and check inequalities in \code{A[11:nrow(A),]} (i.e., all inequalities).
#'
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
#'    \item\code{proportion}: estimated probability that samples are in polytope
#'    \item\code{se}: standard error of probability estimate
#'    \item\code{const_map}: logarithm of the binomial constants that
#'           have to be considered due to equality constraints
#' }
#' @seealso \code{\link{bf_binom}}, \code{\link{count_multinom}}
#'
#' @examples
#' ### a set of linear order constraints:
#' ### x1 < x2 < .... < x6 < .50
#' A <- matrix(c(1, -1, 0, 0, 0, 0,
#'               0, 1, -1, 0, 0, 0,
#'               0, 0, 1, -1, 0, 0,
#'               0, 0, 0, 1, -1, 0,
#'               0, 0, 0, 0, 1, -1,
#'               0, 0, 0, 0, 0, 1),
#'             ncol = 6, byrow = TRUE)
#' b <- c(0, 0, 0, 0, 0, .50)
#'
#' ### check whether a specific vector is inside the polytope:
#' A %*% c(.05, .1, .12, .16, .19, .23) <= b
#'
#'
#' ### observed frequencies and number of observations:
#' k <- c(0, 3, 2, 5, 3, 7)
#' n <- rep(10, 6)
#'
#' ### count prior samples and compare to analytical result
#' prior <- count_binom(0, 0, A, b, M = 5000, steps = 1:4)
#' prior    # to get the proportion: attr(prior, "proportion")
#' (.50)^6 / factorial(6)
#'
#' ### count posterior samples + get Bayes factor
#' posterior <- count_binom(k, n, A, b, M=2000, steps=1:4)
#' count_to_bf(posterior, prior)
#'
#' ### automatic stepwise algorithm
#' prior <- count_binom(0, 0, A, b, M = 500, cmin = 200)
#' posterior <- count_binom(k, n, A, b, M = 500, cmin = 200)
#' count_to_bf(posterior, prior)
#'
#' @template ref_hoijtink2011
#' @template ref_fukuda2004
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
count_binom <- function (k, n, A, b, V, map, prior = c(1, 1), M = 10000,
                         steps, start, cmin = 0, maxiter = 500,
                         burnin = 5, progress = TRUE, cpu = 1){
  check_Mminmax(M, cmin, maxiter, steps)

  if (inherits(cpu, c("SOCKcluster", "cluster")) || is.numeric(cpu) && cpu > 1) {
    arg <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
    count <- run_parallel(arg, fun = "count_binom", cpu = cpu, simplify = "count")
    return(count)
  }

  if (missing(A) || is.null(A)) A <- V
  aggr <- map_k_to_A(k, n, A, map, prior)
  k <- aggr$k
  n <- aggr$n

  if (!missing(b) && !is.null(b)){
    check_Abknprior(A, b, k, n, prior)
    if (cmin == 0 && (missing(steps) || is.null(steps))){
      count <- count_bin(k, n, A, b, prior, M, batch = BATCH, interactive() && progress)

    } else {
      steps <- check_stepsA(steps, A)
      if (length(M) == 2)
        M <- c(M[1], rep(M[2], length(steps) - 1))
      if (missing(start) || is.null(start) || any(start < 0))
        start <-  ml_binom(k + prior[1], n + sum(prior), A, b, map, n.fit = 1, start,
                           control = list(maxit = 50, reltol = .Machine$double.eps^.3))$par
      if (!all(A %*% start < b)){
        # move away into interior
        start <- .5*find_inside(A, b, options = rep(2,length(k)), random = TRUE) + .5*start
      }
      check_start(start, A, b, interior = TRUE)

      if (cmin > 0){
        zeros <- rep(0, length(steps))
        count <- count_auto_bin(k, n, A, b, prior, count = zeros, M = zeros, steps = steps,
                                M_iter = M, cmin = cmin, maxiter = maxiter + length(steps),
                                start, burnin, interactive() && progress)
      } else {
        count <- count_stepwise_bin(k, n, A, b, prior, M, steps, batch = BATCH,
                                    start, burnin, interactive() && progress)
      }
    }
  } else if (!missing(V) && !is.null(V)){
    count <- 0
    m <- M
    a <- c(rbind(k + prior[1], n - k + prior[2]))
    if (interactive() && progress)
      pb <- txtProgressBar(0, M, style = 3)
    while (m > 0 ){
      X <- rpdirichlet(round(BATCH/1000), a, rep(2, ncol(V)), drop_fixed = TRUE)
      count <- count + sum(inside_V(X, V))
      m <- m - round(BATCH/1000)
      if (interactive() && progress)
        setTxtProgressBar(pb, M - m)
    }
    if (interactive() && progress)
      close(pb)
    count <- cbind("count" = count, "M" = M, "steps" = NA)
  } else {
    stop("A/b or V must be provided.")
  }

  attr(count, "const_map_0u") <- aggr$const_map_0u
  as_ineq_count(count)
}

