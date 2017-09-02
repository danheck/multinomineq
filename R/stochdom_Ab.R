#' Ab-Representation for Stochastic Dominance of Histogram Bins
#'
#' Provides the necessary linear equality constraints to test stochastic dominance
#' of continuous distributions, that is, whether the cumulative density functions
#' \eqn{F} satisfy the constraint \eqn{F_1(t) < F_2(t)} for all \eqn{t}.
#'
#' @param bins number of bins of histogram
#' @param conditions number of conditions
#' @param order order constraint on the random variables across conditions.
#'    The default \code{order="<"} implies that the random variables increase across
#'    conditions (implying that the cdfs decrease: \eqn{F_1(t) > F_2(t)}).
#' @seealso \code{\link{stochdom_bf}} to obtain a Bayes factor directly.
#' @template ref_heathcote2010
#' @examples
#' stochdom_Ab(4, 2)
#' stochdom_Ab(4, 3)
#' @export
stochdom_Ab <- function (bins, conditions = 2, order = "<"){
  if (bins != round(bins) | length(bins) != 1 || any(bins < 0))
    stop ("'bins' must be a positive integer.")
  if (conditions != round(conditions) | length(conditions) != 1 || any(conditions < 2))
    stop ("'conditions' must be an integer larger than 1.")

  dnames <- c(t(outer(paste0("c", 1:conditions),
                      paste0("b", 1:(bins-1)),
                      paste, sep = ",")))
  A <- matrix(0, 0, (bins - 1) * conditions,
              dimnames = list(NULL, dnames))
  zeros <- plus <- minus <- matrix(0, bins-1, bins-1)
  plus[lower.tri(zeros, diag = TRUE)]  <- +1
  minus[lower.tri(zeros, diag = TRUE)] <- -1

  for (c in 1:(conditions - 1)){
    D <- matrix(0, bins - 1, ncol(A),
                dimnames = list(paste0("c",c,order,"c",c+1, ",",1:(bins-1)), NULL))
    idx1 <- (c-1)*(bins-1) + 1:(bins-1)
    idx2 <-     c*(bins-1) + 1:(bins-1)
    if (order == ">"){
      D[,idx1] <- plus
      D[,idx2] <- minus
    } else if (order == "<") {
      D[,idx1] <- minus
      D[,idx2] <- plus
    }
    A <- rbind(A, D)
  }
  list("A" = A, "b" = rep(0, nrow(A)))
}
