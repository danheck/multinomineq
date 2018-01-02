#' Bayes Factor for Nonlinear Order Constraints
#'
#' Draws prior and posterior samples from product-multinomial model and computes the
#' Bayes factor for a user-specified nonlinear order constraint defined via an
#' indicator function of the parameters: (p11,p12,p13,  p21,p22,...)
#'
#' @inheritParams count_multinom
#' @inheritParams count_binom
#' @param inside a function that takes a vector with probabilities
#'    \code{p=c(p11,p12,p13,  p21,p22,...)} as input and returns \code{1} or \code{TRUE}
#'    if the order constraints are satisfied and \code{0} or \code{FALSE} if not.
#' @template ref_klugkist2007
#' @examples
#' ### 2x2x2 continceny table (Klugkist & Hojtink, 2007)
#' # (defendant's race) x (victim's race) x (death penalty)
#' # index: 0 = white/white/yes  ; 1 = black/black/no
#'
#' # p = (p000,p001,  p010,p011,  p100,p101,  p110,p111)
#' k <- c(19,132,     0, 9,       11, 52,     6,97)
#'
#' # Model2: p000*p101 < p100*p001 & p010*p111 < p110*p011
#' model2 <- function(x)
#'    x[1]*x[6] < x[5]*x[2] & x[3]*x[8] < x[7]*x[4]
#' bf_nonlinear(k, 8, model2, M=1e5)  # K&H07: bf_0u=1.62
#' @export
bf_nonlinear <- function(k, options, inside,
                         prior = rep(1, sum(options)), M = 10000){
  check_io(inside, options)
  check_Mminmax(M)
  todo <- M
  cnt_po <- cnt_pr <- 0
  while (todo > 0){
    po <- rpdirichlet(min(BATCH, todo), k + prior, options, p_drop = FALSE)
    pr <- rpdirichlet(min(BATCH, todo), prior, options, p_drop = FALSE)
    cnt_po <- cnt_po + sum(apply(po, 1, inside))
    cnt_pr <- cnt_pr + sum(apply(pr, 1, inside))
    todo <- todo - BATCH
  }
  count_to_bf(cbind("count" = cnt_po, "M" = M),
              cbind("count" = cnt_pr, "M" = M))
}
