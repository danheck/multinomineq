#' Automatic Construction of Ab-Representation for Common Inequality Constraints
#'
#' Constructs the matrix \code{A} and vector \code{b} of the Ab-representation
#' \code{A*x < b} for common inequality constraints such as "the probability j is
#' larger than all others (\code{Ab_max})" or "the probabilities are ordered
#' (\code{Ab_monotonicity})").
#'
#' @param which_max vector of indices refering to probabilities that are
#'   assumed to be larger than the remaining probabilities
#'   (e.g., \code{which_max=c(1,2)} means that \code{p1>p3, p1>p4,..., p2>p3, ...}).
#'   Note that the indices refer to \emph{all} probabilities/categories (including
#'   one fixed probability within each multinomial distribution).
#' @param exclude vector of indices refering to probabilities that are
#'     excluded from the construction of the order constraints (including
#'     the fixed probabilities).
#' @param exclude_fixed whether to exclude the fixed probabilities (i.e., the
#'   last probability within each multinomial) from the construction of the
#'   order constraints. For example, if \code{options=c(2,2,3)} then the
#'   probabilities/columns 2, 4, and 7 are dropped (which is equivalent to
#'   \code{exclude=c(2,4,7)}). This option is usually appropriate for binomial
#'   probabilities (i.e., if \code{options = c(2,2,2,...)}), e.g., when the
#'   interest is in the probability of correct responding across different item types.
#' @param drop_fixed whether to drop columns of \code{A} containing the
#'   fixed probabilities (i.e., the last probability within each multinomial).
#'   \emph{after} construction of the inequalities.
#' @inheritParams count_multinom
#'
#' @return a list with the matrix \code{A} and the vectors \code{b} and \code{options}
#'
#' @examples
#' # Example 1: Multinomial with 5 categories
#' # Hypothesis: p1 is larger than p2,p3,p4,p5
#' Ab_max(which_max = 1, options = 5)
#'
#' # Example 2: Four binomial probabilities
#' # Hypothesis: p1 is larger than p2,p3,p4
#' Ab_max(which_max = 1, options = c(2,2,2,2), exclude_fixed = TRUE)
#' @export
Ab_max <- function(which_max, options, exclude = c(), exclude_fixed = FALSE,
                   drop_fixed = TRUE){

  check_o(options)
  stopifnot(length(which_max) < sum(options) - 1,
            round(which_max) == which_max)
  stopifnot(length(exclude) <= sum(options) - length(which_max) - 1,
            is.null(exclude) || all(round(exclude) == exclude))

  if (exclude_fixed)
    exclude <- sort(unique(c(exclude, cumsum(options))))

  J <- sum(options)
  J_max <- length(which_max)
  J_exc <- length(exclude)
  smaller <- setdiff(1:J, c(which_max, exclude))

  rows_per_max <- J - J_max - length(exclude)
  A <- matrix(0, rows_per_max * length(which_max), J)
  colnames(A) <- index_mult(options, fixed = TRUE)
  rownames(A) <- 1:nrow(A)
  b <- rep(0, nrow(A))

  for (mm in seq(which_max)){
    row_idx <- (mm - 1)*rows_per_max + 1:rows_per_max
    A[row_idx, which_max[mm]] <- -1
    A[row_idx, smaller] <- diag(length(smaller))
    rownames(A)[row_idx] <- paste0(colnames(A)[smaller], "<",
                                   colnames(A)[which_max[mm]])
  }

  Ab <- list(A = A, b = b, options = options)
  if (drop_fixed)
    Ab <- Ab_drop_fixed(A, b, options)
  Ab
}
