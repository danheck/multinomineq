#' Posterior Sampling for Inequality-Constrained Multinomial Models
#'
#' Uses Gibbs sampling to draw posterior samples for binomial and multinomial
#' models with linear inequality-constraints.
#'
#' @inheritParams inside
#' @inheritParams count_multinom
#' @inheritParams count_binom
#' @inheritParams sampling_binom
#' @param M number of posterior samples
#' @param burnin number of burnin samples that are discarded. Can be chosen to be
#'     small if the maxmimum-likelihood estimate is used as a starting value.
#' @param cpu either the number of CPUs using separate MCMC chains in parallel,
#'     or a parallel cluster (e.g., \code{cl <- parallel::makeCluster(3)}).
#'     All arguments of the function call are passed directly to each core,
#'     and thus the total number of samples is \code{M*number_cpu}.
#'
#' @details
#' Draws posterior samples for binomial/multinomial random utility models that
#' assume a mixture over predefined preference orders/vertices that jointly define
#' a convex polytope via the set of inequalities \code{A * x <= b} or as the
#' convex hull of a set of vertices \code{V}.
#'
#' @seealso \code{\link{count_multinom}}, \code{\link{ml_multinom}}
#' @return an \code{mcmc} matrix (or an \code{mcmc.list} if \code{cpu>1}) with
#'     posterior samples of the binomial/multinomial probability parameters.
#'     See \code{\link[coda]{mcmc}}) .
#' @template ref_myung2005
#'
#' @examples
#' ############### binomial ##########################
#' A <- matrix(c(1, 0, 0,   # x1 < .50
#'               1, 1, 1,   # x1+x2+x3 < 1
#'               0, 2, 2,   # 2*x2+2*x3 < 1
#'               0, -1, 0,  # x2 > .2
#'               0, 0, 1),  # x3 < .1
#'             ncol = 3, byrow = TRUE)
#' b <- c(.5, 1, 1, -.2, .1)
#' samp <- sampling_binom(c(5,12,7), c(20,20,20), A, b)
#' head(samp)
#' apply(samp, 2, plot, type = "l", ylim = c(0,.6))
#'
#'
#' ############### multinomial ##########################
#' # binary and ternary choice:
#' #           (a1,a2   b1,b2,b3)
#' k       <- c(15,9,   5,2,17)
#' options <- c(2,      3)
#'
#' # columns:   (a1,  b1,b2)
#' A <- matrix(c(1, 0, 0,   # a1 < .20
#'               0, 2, 1,   # 2*b1+b2 < 1
#'               0, -1, 0,  # b1 > .2
#'               0, 0, 1),  # b2 < .4
#'             ncol = 3, byrow = TRUE)
#' b <- c(.2, 1, -.2, .4)
#' samp <- sampling_multinom(k, options, A, b)
#' head(samp)
#' apply(samp, 2, plot, type = "l", ylim = c(0, 1))
#' @importFrom coda as.mcmc as.mcmc.list mcmc mcmc.list
#' @export
sampling_multinom <- function (k, options, A, b, V, prior = rep(1, sum(options)),
                               M = 5000, start, burnin = 10,
                               progress = TRUE, cpu = 1){
  if (missing(k) || is.null(k) || (length(k) == 1 && k == 0))
    k <- rep(0, sum(options))

  if (class(cpu) %in% c("SOCKcluster", "cluster") || is.numeric(cpu) && cpu > 1) {
    arg <- lapply(as.list(match.call())[-1], eval, envir = environment())
    mcmc_list <- run_parallel(arg, fun = "sampling_multinom", cpu = cpu, simplify = "as.mcmc.list")
    return(mcmc_list)
  }

  if (!missing(V) && !is.null(V)){
    mcmc <- sampling_V(k, options, V, prior = prior, M = M,
                       start = start, burnin = burnin, progress = progress)
  } else {
    if (missing(start) || is.null(start) || any(start < 0))
      start <-  ml_multinom(k, options, A, b, n.fit = 1, maxit = 20)$par
    check_start(start, A, b, interior = TRUE)
    if (length(prior) == 1) prior <- rep(prior, sum(options))
    check_Abokprior(A, b, options, k, prior)
    mcmc <- sampling_mult(k, options, A, b, prior, M, start, burnin, progress)
    colnames(mcmc) <- colnames(A)
    if (is.null(colnames(A)))
      colnames(mcmc) <- drop_fixed(index_mult(options), options)
  }
  mcmc(mcmc, start = burnin + 1, end = M)
}


#' @rdname sampling_multinom
#' @export
sampling_binom <- function (k, n, A, b, V, map = 1:ncol(A), prior = c(1, 1),
                            M = 5000, start, burnin = 10,
                            progress = TRUE, cpu = 1){
  if (length(n) == 1) n <- rep(n, length(k))
  if (length(prior) == 1) prior <- rep(prior, length(k))

  if (class(cpu) %in% c("SOCKcluster", "cluster") || is.numeric(cpu) && cpu > 1) {
    arg <- lapply(as.list(match.call())[-1],
                  function(i) tryCatch(eval(i), error = function(e) NULL))
    mcmc.list <- run_parallel(arg, fun = "sampling_binom", cpu = cpu, simplify = "as.mcmc.list")
    return(mcmc.list)
  }

  if (!missing(V) && !is.null(V)){
    options <- rep(2, length(k))
    k2 <- add_fixed(k, options = options, sum = n)
    mcmc <- sampling_V(k2, options = options, V = V,
                       prior = rep(prior, length(k)), M = M,
                       start = start, burnin = burnin, progress = progress)
  } else {
    aggr <- map_k_to_A(k, n, A, map)
    k <- aggr$k
    n <- aggr$n
    check_Abknprior(A, b, k, n, prior)
    if (missing(start) || is.null(start) || any(start < 0))
      start <-  ml_binom(k, n, A, b, map, n.fit = 1, maxit = 20)$par
    check_start(start, A, b, interior = TRUE)
    mcmc <- sampling_bin(k, n, A, b, prior, M, start, burnin, progress)
    colnames(mcmc) <- colnames(A)
  }
  mcmc(mcmc, start = burnin + 1, end = M)
}
