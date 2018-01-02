#' Aggregation of Individual Bayes Factors
#'
#' Aggregation of multiple individual (N=1) Bayes factors to obtain the evidence
#' for a hypothesis in a population of persons.
#'
#'
#' @param bfs a vector with individual Bayes factors,
#'     a matrix with one type of Bayes-factor comparison per column,
#'     or a list of matrices with a named column \code{"bf"} (as returned by
#'     \code{\link{bf_multinom}}/\code{\link{count_to_bf}}).
#' @return a vector or matrix with named elements/columns:
#'    \itemize{
#'      \item \code{population_bf}: the product of individual BFs
#'      \item \code{geometric_bf}: the geometric mean of the individual BFs
#'      \item \code{evidence_rate}: the proportion of BFs>1 (BFs<1) if \code{geometric_bf>1} (<1).
#'          Values close to 1.00 indicate homogeneity.
#'      \item \code{stability_rate}: the proportion \code{bfs>geometric_bf} (<) if \code{geometric_bf>1} (<).
#'          Values close to 0.50 indicate stability.
#'    }
#' @template ref_klaassen2018
#'
#' @examples
#' # consistent evidence across persons:
#' bfs <- c(2.3, 1.8, 3.3, 2.8, 4.0, 1.9, 2.5)
#' population_bf(bfs)
#'
#' # (A) heterogeneous, inconsistent evidence
#' # (B) heterogeneous, inconsistent evidence
#' bfs <- cbind(A = c(2.3, 1.8, 3.3, 2.8, 4.0, 1.9, 2.5),
#'              B = c(10.3, .7, 3.3, .8, 14.0, .9, 1.5))
#' population_bf(bfs)
#'
#' @export
population_bf <- function(bfs){
  UseMethod("population_bf", bfs)
}

#' @export
population_bf.list <- function(bfs){
  check_bf <-
    all(sapply(bfs, is.matrix)) &&
    all("bf" %in% sapply(bfs, colnames)) &&
    all(sapply(bfs, nrow) == nrow(bfs[[1]]))  &&
    all(sapply(bfs, rownames) == rownames(bfs[[1]]))

  if (!check_bf)
    stop("If a list is supplied to 'population_bf', each element must be a matrix with\n:",
         "   - a column named 'bf' (See ?count_to_bf)\n",
         "    - identical number of rows and rownames (same BFs per person)")
  bf_mat <- t(sapply(bfs, function(b) b[,"bf"]))
  population_bf(bf_mat)
}

#' @export
population_bf.matrix <- function(bfs){
  apply(bfs, 2, population_bf)
}

#' @export
population_bf.data.frame <- function(bfs){
  apply(bfs, 2, population_bf)
}

#' @export
population_bf.numeric <- function(bfs){
  stopifnot(all(bfs >= 0))
  p_bf <- prod(bfs)
  gp_bf <- p_bf^(1 / length(bfs))
  # evidence rate
  er <- ifelse(gp_bf < 1, mean(bfs < 1), mean(bfs > 1))
  # stability rate
  sr <- ifelse(gp_bf < 1, mean(bfs < gp_bf), mean(bfs > gp_bf))

  c(population_bf = p_bf, geometric_bf = gp_bf,
    evidence_rate = er, stability_rate = sr)
}

#' @export
population_bf.default <- function(bfs){
  stop("'bfs' must be numeric.")
}
