#' Bayes Factor for Nonlinear Inequality Constraints
#'
#' Computes the encompassing Bayes factor for a user-specified, nonlinear inequality
#' constraint. Restrictions are defined via an indicator function of the free parameters
#' \code{c(p11,p12,p13,  p21,p22,...)} (i.e., the multinomial probabilities).
#'
#' @inheritParams bf_multinom
#' @inheritParams count_multinom
#' @inheritParams count_binom
#' @param inside an indicator function that takes a vector with probabilities
#'    \code{p=c(p11,p12,  p21,p22,...)} (where the last probability for each
#'    multinomial is dropped) as input and returns \code{1} or \code{TRUE}
#'    if the order constraints are satisfied and \code{0} or \code{FALSE} otherwise
#'    (see details).
#'
#' @template details_indicator
#' @template ref_klugkist2007
#' @template ref_klugkist2010
#'
#' @examples
#' ##### 2x2x2 continceny table (Klugkist & Hojtink, 2007)
#' #
#' # (defendant's race) x (victim's race) x (death penalty)
#' # indexing: 0 = white/white/yes  ; 1 = black/black/no
#' # probabilities: (p000,p001,  p010,p011,  p100,p101,  p110,p111)
#' # Model2:
#' # p000*p101 < p100*p001  &   p010*p111 < p110*p011
#'
#' # observed frequencies:
#' k <- c(19, 132, 0, 9, 11, 52, 6, 97)
#'
#' model <- function(x) {
#'   x[1] * x[6] < x[5] * x[2] & x[3] * (1 - sum(x)) < x[7] * x[4]
#' }
#' # NOTE: "1-sum(x)"  must be used instead of "x[8]"!
#'
#' # compute Bayes factor (Klugkist 2007: bf_0u=1.62)
#' bf_nonlinear(k, 8, model, M = 50000)
#'
#' \donttest{
#' ##### Using a C++ indicator function (much faster)
#' cpp_code <- "SEXP model(NumericVector x){
#'   return wrap(x[0]*x[5] < x[4]*x[1] & x[2]*(1-sum(x)) < x[6]*x[3]);}"
#' # NOTE: C++ indexing starts at 0!
#'
#' # define C++ pointer to indicator function:
#' model_cpp <- RcppXPtrUtils::cppXPtr(cpp_code)
#'
#' bf_nonlinear(
#'   k = c(19, 132, 0, 9, 11, 52, 6, 97), M = 100000,
#'   options = 8, inside = model_cpp
#' )
#' }
#' @export
bf_nonlinear <- function(k, options, inside, prior = rep(1, sum(options)),
                         log = FALSE, ...) {
  arg <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  arg$log <- NULL
  po <- do.call("count_nonlinear", arg)
  arg$k <- 0
  pr <- do.call("count_nonlinear", arg)
  count_to_bf(po, pr, log = log)
}



#' @rdname bf_nonlinear
#' @export
count_nonlinear <- function(k = 0, options, inside, prior = rep(1, sum(options)),
                            M = 5000, progress = TRUE, cpu = 1) {
  check_io(inside, options)
  stopifnot(length(M) == 1, M > 0, M == round(M))

  if (inherits(cpu, c("SOCKcluster", "cluster")) || is.numeric(cpu) && cpu > 1) {
    if (inherits(inside, "XPtr")) {
      stop("C++ functions (defined via RcppXPtrUtils::cppXPtr) not supported if cpu>1.")
    }
    arg <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
    count <- run_parallel(arg, fun = "count_nonlinear", cpu = cpu, simplify = "count")
    return(count)
  }

  if (length(k) == 1 && k == 0) k <- rep(0, sum(options))
  if (length(prior) == 1) prior <- rep(prior, sum(options))

  if (is.function(inside)) {
    todo <- M
    cnt <- 0
    while (todo > 0) {
      pp <- rpdirichlet(min(BATCH, todo), k + prior, options)
      count <- cnt + sum(apply(pp, 1, inside))
      todo <- todo - BATCH
    }
  } else {
    count <- count_nonlin_cpp(
      k, options, inside, prior, M, BATCH,
      interactive() && progress
    )
  }

  as_ineq_count(cbind(count = count, M = M, steps = 1))
}
