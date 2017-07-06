#' Sort Inequalities by Predictive Strength
#'
#' Uses uniform sampling to order the inequalities by their predictive strength.
#' The constraint that constrains the paramters most is placed at the first position.
#'
#' @inheritParams compute_bf
#' @param M number of uniform samples from the unit cube
#' @examples
#' b <- c(0,0,.30,.70)
#'  A <- matrix(c(-1,1,0,  # p1 >= p2
#'                0,1,-1,  # p2 <= p3
#'                1,0,0,   # p1 <=.30
#'                0,1,0),  # p2 <= .70
#'                ncol = 3, byrow = 2)
#' sort_inequalities(A, b)
#' @export
sort_inequalities <- function (A, b, M = 1e5){
  check_kAb(rep(1, ncol(A)), A, b)
  S <- ncol(A)
  x <- matrix(runif(M * S), S)
  test <- A %*% x <= b

  volume = rowMeans(test)
  o <- order(volume)
  list("A" = A[o,], "b" = b[o], "volume" = volume[o])
}
