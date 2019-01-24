#' Count How Many Samples Satisfy Linear Inequalities (Multinomial)
#'
#' Draws prior/posterior samples for product-multinomial data and counts how many samples are
#' inside the convex polytope defined by
#' (1) the inequalities \code{A*x <= b} or
#' (2) the convex hull over the vertices V.
#'
#' @param A a matrix defining the convex polytope via \code{A*x <= b}.
#'    The columns of \code{A} do not include the last choice option per item type and
#'    thus the number of columns must be equal to \code{sum(options-1)}
#'    (e.g., the column order of \code{A} for \code{k = c(a1,a2,a2, b1,b2)}
#'    is \code{c(a1,a2, b1)}).
#' @param options the number of observable categories per item type,
#'    e.g., \code{c(3,2)} for a ternary and binary item.
#'     The sum of \code{options} must be equal to the length of \code{k}.
#' @param k the number of choices for each alternative ordered by item type (e.g.
#'     \code{c(a1,a2,a3,  b1,b2)} for a ternary and a binary item type).
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
#'    \item\code{proportion}: estimated probability that samples are in polytope
#'    \item\code{se}: standard error of probability estimate
#' }
#' @template ref_hoijtink2011
#' @seealso \code{\link{bf_multinom}}, \code{\link{count_binom}}
#'
#' @examples
#' ### frequencies:
#' #           (a1,a2,a3, b1,b2,b3,b4)
#' k <-       c(1,5,9,    5,3,7,8)
#' options <- c(3,        4)
#'
#' ### linear order constraints
#' # a1<a2<a3   AND   b2<b3<.50
#' # (note: a2<a3 <=> a2 < 1-a1-a2 <=> a1+2*a2 < 1)
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
#' bf_multinom(k, options, A, b, M=10000)
#' bf_multinom(k, options, A, b, cmin = 5000, M = 1000)
#' @export
count_multinom <- function (k = 0, options, A, b, V, prior = rep(1, sum(options)),
                            M = 5000, steps, start, cmin = 0, maxiter = 500,
                            burnin = 5, progress = TRUE, cpu = 1){

  if (class(cpu) %in% c("SOCKcluster", "cluster") || is.numeric(cpu) && cpu > 1) {
    arg <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
    count <- run_parallel(arg, fun = "count_multinom", cpu = cpu, simplify = "count")
    return(count)
  }

  if (length(k) == 1 && k == 0) k <- rep(0, sum(options))
  if (length(prior) == 1) prior <- rep(prior, sum(options))

  if (!missing(b) && !is.null(b)){
    check_Abokprior(A, b, options, k, prior)
    check_Mminmax(M, cmin, maxiter)

    if (cmin == 0 && (missing(steps) || is.null(steps))){
      count <- count_mult(k, options, A, b, prior, M, batch = BATCH, progress)

    } else {
      steps <- check_stepsA(steps, A)
      if (missing(start) || is.null(start) || any(start < 0))
        start <-  ml_multinom(k + prior, options, A, b, n.fit = 1, start,
                              control = list(maxit = 50, reltol = .Machine$double.eps^.3))$par
      check_start(start, A, b, interior = FALSE)

      if (cmin > 0){
        zeros <- rep(0, length(steps))
        count <- count_auto_mult(k, options, A, b, prior, zeros, zeros, steps, ## SUM TO ZERO!!!!
                                 M_iter = M, cmin = cmin, maxiter = maxiter + length(steps),
                                 start, burnin, progress)
      } else {
        count <- count_stepwise_multi(k, options, A, b, prior, M, steps,
                                      batch = BATCH, start, burnin, progress)
      }
    }

  } else if (!missing(V) && !is.null(V)){
    check_Vx(V, drop_fixed(k, options))
    count <- 0
    m <- M
    a <- k + prior #c(rbind(k + prior[1], n - k + prior[2]))
    if (progress) pb <- txtProgressBar(0, M, style = 3)
    while (m > 0 ){
      X <- rpdirichlet(n = round(BATCH/1000), alpha = a, options = options, p_drop = TRUE)
      count <- count + sum(inside_V(X, V))
      m <- m - round(BATCH/1000)
      if (progress) setTxtProgressBar(pb, M - m)
    }
    if (progress) close(pb)
    count <- cbind("count" = count, "M" = M, "steps" = NA)

  } else {
    stop("A/b or V must be provided.")
  }
  as_ineq_count(count)
}
