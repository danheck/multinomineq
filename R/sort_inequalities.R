#' Sort Inequalities by Acceptance Rate
#'
#' Uses samples from the prior/posterior to order the inequalities by the acceptance rate.
#'
#' @inheritParams inside
#' @param k optional: number of observed frequencies (only for posterior sampling).
#' @param options optional: number of options per item type/category system.
#'      Uniform sampling on [0,1] for each parameter is used if omitted.
#' @param M number of samples.
#' @param drop_irrelevant whether to drop irrelevant constraints for probabilities such as
#'     \code{theta[1] >= 0}, \code{theta[1] <= 1}, or \code{sum(theta) <= 1}.
#' @details
#'
#' Those constraints that are rejected most often are placed at the first positions.
#' This can help when computing the encompassing Bayes factor and counting how many samples
#' satisfy the constraints (e.g., \code{\link{count_binom}} or \code{\link{bf_multinom}}).
#' Essentially, it becomes more likely that the while-loop for testing
#' whether the inequalities hold can stop earlier, thus making the computation faster.
#'
#' The function could also be helpful to improve the efficiency of the stepwise
#' sampling implemented in \code{\link{count_binom}} and \code{\link{count_multinom}}.
#' First, one can use accept-reject sampling to test the first few, rejected
#' inequalities. Next, one can use a Gibbs sampler to draw samples conditional on the
#' first constraints.
#'
#'
#' @examples
#' ### Binomial probabilities
#' b <- c(0,0,.30,.70, 1)
#' A <- matrix(c(-1,1,0,  # p1 >= p2
#'               0,1,-1,  # p2 <= p3
#'               1,0,0,   # p1 <=.30
#'               0,1,0,   # p2 <= .70
#'               0,0,1),  # p3 <= 1 (redundant)
#'               ncol = 3, byrow = 2)
#' Ab_sort(A, b)
#'
#'
#' ### Multinomial probabilities
#' # prior sampling:
#' Ab_sort(A, b, options = 4)
#' # posterior sampling:
#' Ab_sort(A, b, k = c(10,3, 2, 14), options = 4)
#'
#' @export
Ab_sort <- function (A, b, k = 0, options, M = 1000, drop_irrelevant = TRUE){
  check_Ab(A, b)
  S <- ncol(A)
  if (missing(options)){
    x <- matrix(runif(M * S), M, S)
  } else {
    if (length(k) == 1)
      k <- rep(k, sum(options))
    check_ko(k, options)
    x <- rpdirichlet(M, k + 1, options)
  }

  if (drop_irrelevant){
    Ab <- Ab_drop_irrelevant(A, b, options)
    A <- Ab$A
    b <- Ab$b
  }
  accept <- A %*% t(x) <= b
  accept_rate <- rowMeans(accept)
  o <- order(accept_rate, decreasing = FALSE)

  list("A" = A[o,], "b" = b[o], "accept_rate" = accept_rate[o])
}


### drop constraints: 0<p<1
Ab_drop_irrelevant <- function(A, b, options){
  lower1   <- apply(A, 1, function(a) sum(a==0) == ncol(A)-1 && sum(a==1) == 1) & b==1
  greater0 <- apply(A, 1, function(a) sum(a==0) == ncol(A)-1 && sum(a==-1) == 1) & b==0
  sum1 <- rep(FALSE, nrow(A))
  if (!missing(options)){
    for (i in seq_along(options)){
      prev <- sum(options[seq(0, i-1)] - 1)
      idx <- prev + seq(options[i] - 1)
      checki <- b == 1 & apply(A[,idx,drop=FALSE]==1, 1, all) & apply(A[,-idx,drop=FALSE]==0, 1, all)
      sum1 <- sum1 | checki
    }
  }
  A1 <- A[!lower1 & !greater0 & !sum1,,drop = FALSE]
  b1 <- b[!lower1 & !greater0 & !sum1]

  list(A = A1, b = b1)
}

# ### split A into block-diagonal matrices => allows to compute encompassing BF independently
# Ab_split <- function(A, b){
#
#   order1 <- order(-abs(A1[,1]))
#   A2 <- A1[order1,]
#   b2 <- b1[order1]
#   for (i in 2:ncol(A2) ){
#     sel <- apply(A2[,seq(1, i - 1), drop = FALSE] == 0, 1, all)
#     if (any(sel)){
#       orderi <- order(-abs(A2[sel,i]))
#       A2[sel,] <- A2[sel,,drop=FALSE][orderi,]
#       b2[sel] <- b2[sel][orderi]
#     }
#   }
#   list(cbind(A2, .....b2 = b2))
# }
