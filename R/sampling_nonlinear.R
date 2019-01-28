#' Posterior Sampling for Multinomial Models with Nonlinear Inequalities
#'
#' A Gibbs sampler that draws posterior samples of probability parameters
#' conditional on a (possibly nonlinear) indicator function defining a
#' restricted parameter space that is convex.
#'
#' @param eps precision of the bisection algorithm
#' @inheritParams bf_nonlinear
#' @inheritParams sampling_multinom
#'
#' @template details_indicator
#'
#' @details
#' For each parameter, the Gibbs sampler draws a sample from the
#' conditional posterior distribution (a scaled, truncated beta).
#' The conditional truncation boundaries are computed with a bisection algorithm.
#' This requires that the restricted parameteter space defined by the indicator
#' function is convex.
#'
#' @examples
#' # two binomial success probabilities: x = c(x1, x2)
#' # restriction to a circle:
#' model <- function(x)
#'   (x[1]-.50)^2 + (x[2]-.50)^2 <= .15
#'
#' # draw prior samples
#' mcmc <- sampling_nonlinear(k = 0, options = c(2,2),
#'                            inside = model, M = 1000)
#' head(mcmc)
#' plot(c(mcmc[,1]), c(mcmc[,2]), xlim=0:1, ylim=0:1)
#'
#'
#'
#' \dontrun{
#' ##### Using a C++ indicator function (much faster)
#' cpp_code <- "SEXP inside(NumericVector x){
#'   return wrap( sum(pow(x-.50, 2)) <= .15);}"
#' # NOTE: Uses Rcpp sugar syntax (vectorized sum & pow)
#'
#' # define function via C++ pointer:
#' model_cpp <- RcppXPtrUtils::cppXPtr(cpp_code)
#' mcmc <- sampling_nonlinear(k=0, options = c(2,2),
#'                            inside = model_cpp, M=1000)
#' head(mcmc)
#' plot(c(mcmc[,1]), c(mcmc[,2]), xlim=0:1, ylim=0:1)
#' }
#' @export
sampling_nonlinear <- function (k, options, inside, prior = rep(1, sum(options)),
                                M = 1000, start, burnin = 10, eps = 1e-6,
                                progress = TRUE, cpu = 1){

  stopifnot(M > burnin, burnin > 0)
  if (missing(k) || is.null(k) || (length(k) == 1 && k == 0))
    k <- rep(0, sum(options))
  check_ko(k, options)
  check_io(inside, options)

  if (class(cpu) %in% c("SOCKcluster", "cluster") || is.numeric(cpu) && cpu > 1) {
    if (class(inside) == "XPtr")
      stop("C++ functions (defined via RcppXPtrUtils::cppXPtr) not supported if cpu>1.")
    arg <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
    if (is.null(arg$k) || (length(arg$k) == 1 && arg$k == 0))
      arg$k <- rep(0, sum(arg$options))
    mcmc_list <- run_parallel(arg, fun = "sampling_nonlinear", cpu = cpu, simplify = "as.mcmc.list")
    return(mcmc_list)
  }
  if (length(prior) == 1) prior <- rep(prior, sum(options))

  # find random starting values (and use the one with highest posterior density)
  if (missing(start)  || is.null(start) || anyNA(start)){
    start <- NA
    cnt <- 0
    while (anyNA(start) && cnt < 100){
      proposal <- rpdirichlet(100, rep(1, sum(options)), options)
      if (is.function(inside))
        idx <- which(apply(proposal, 1, inside))
      else
        idx <- which(apply(proposal, 1,
                           function(x) call_xptr(inside, x)) == 1)
      if (length(idx) > 0){
        ll <- loglik_multinom(proposal[idx,,drop = FALSE], k + prior, options)
        start <- proposal[idx[which.min(ll)],]

        # rough convexity check:
        alpha <- rpdirichlet(100, rep(1, length(idx)), length(idx), p_drop = FALSE)
        xx <- alpha %*% proposal[idx,,drop=FALSE]
        if (is.function(inside))
          convex <- apply(xx, 1, inside)
        else
          convex <- apply(xx, 1, function(x) call_xptr(inside, x))
        if (!all(convex == 1))
          stop("Indicator function 'inside' does not define a convex parameter space!\n",
               "  (note that multiplicative constraints such as x[1]*x[2]<0.50 are not convex)")
      }
    }
    if (anyNA(start))
      stop("Could not find starting values (based on 10000 samples).")
  }
  stopifnot(length(start) == sum(options - 1))
  stopifnot(all(tapply(start, rep(seq(options), options-1), sum) <= 1))

  if (is.function(inside)){
    stopifnot(inside(start))
    mcmc <- sampling_nonlin_r(k, options, inside, prior, M, start, burnin, progress, eps)
  } else {
    stopifnot(call_xptr(inside, start) == 1)
    mcmc <- sampling_nonlin_cpp(k, options, inside, prior, M, start, burnin, progress, eps)
  }

  colnames(mcmc) <- index_mult(options, fixed = FALSE)
  mcmc(mcmc, start = burnin + 1, end = M)
}





