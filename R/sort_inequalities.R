#' Sort Inequalities by Predictive Strength
#'
#' Uses uniform sampling to order the inequalities by their predictive strength.
#' The constraint that restricts the paramter space most strongly is placed at the first position.
#'
#' @inheritParams inside
#' @param M number of uniform samples from the unit cube
#' @details This function might be helpful to improve the efficiency of the stepwise
#'     sampling implemented in \code{\link{count_binomial}} and \code{\link{count_multinomial}}.
#'
#' @examples
#' b <- c(0,0,.30,.70)
#'  A <- matrix(c(-1,1,0,  # p1 >= p2
#'                0,1,-1,  # p2 <= p3
#'                1,0,0,   # p1 <=.30
#'                0,1,0),  # p2 <= .70
#'                ncol = 3, byrow = 2)
#' sort_Ab(A, b)
#' @export
sort_Ab <- function (A, b, M = 1000){
  check_Ab(A, b)
  S <- ncol(A)
  x <- matrix(runif(M * S), S)
  test <- A %*% x <= b

  volume = rowMeans(test)
  o <- order(volume, decreasing = FALSE)
  list("A" = A[o,], "b" = b[o], "volume" = volume[o])
}
