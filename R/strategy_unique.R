#' Unique Patterns/Item Types of Strategy Predictions
#'
#' Find unique item types, which are defined as patterns of cue values
#' that lead to identical strategy predictions.
#'
#' @param strategies a list of strategy predictions with the same length of
#'     the vector \code{pattern}, see \link{strategy_multiattribute}.
#' @param add_baseline whether to add a baseline model which assumes one probability in [0,1] for each item type.
#' @param reversed whether reversed patterns are treated separately
#'    (default: automatically switch Option A and B if \code{pattern=c(-1,1,1,1)})
#' @return a list including:
#' \itemize{
#'    \item \code{unique}: a matrix with the unique strategy patterns
#'    \item \code{item_type}: a vector that maps the original predictions to item types (negative: reversed items)
#'    \item \code{strategies}: a list with strategy predictions with \code{pattern} adapted to the unique item types
#'  }
#'
#' @examples
#' data(heck2017_raw)
#' ca <- heck2017_raw[1:100, c("a1", "a2", "a3", "a4")]
#' cb <- heck2017_raw[1:100, c("b1", "b2", "b3", "b4")]
#' v <- c(.9, .8, .7, .6)
#' strats <- strategy_multiattribute(
#'   ca, cb, v,
#'   c("WADDprob", "WADD", "TTB")
#' )
#' strategy_unique(strats)
#' @export
strategy_unique <- function(strategies,
                            add_baseline = TRUE, reversed = FALSE) {
  sapply(strategies, check_strategy)
  n_strat <- length(strategies)
  if (is.null(names(strategies))) {
    for (s in 1:n_strat) {
      if (!is.null(strategies[[s]]$label)) {
        names(strategies)[s] <- strategies[[s]]$label
      }
    }
  }
  patterns <- matrix(sapply(strategies, "[[", "pattern"), ncol = n_strat)
  colnames(patterns) <- names(strategies)
  p_sign <- sign(patterns)
  n_unique <- apply(patterns, 2, function(p) length(unique))
  if (max(n_unique) == nrow(patterns)) {
    warning("Number of item types cannot be reduced. Check that 'baseline' model is not included in list.")
  }

  p_unique <- unique(patterns)
  if (!reversed) {
    p_rev <- -p_unique
    p_reduced <- p_unique[p_unique[, 1] >= 0, ]
    idx_rev <- match_pattern(-p_reduced, p_unique)
    p_unique <- p_reduced
  }
  rownames(p_unique) <- paste0("item_type_", 1:nrow(p_unique))
  item_type <- match_pattern(patterns, p_unique)
  if (!reversed) {
    sel_rev <- is.na(item_type)
    item_type[sel_rev] <- -match_pattern(-patterns[sel_rev, ], p_unique)
  }
  for (s in 1:n_strat) {
    strategies[[s]]$pattern <- p_unique[, s]
  }
  if (add_baseline) {
    strategies$baseline <- as_strategy(1:nrow(p_unique), 1, FALSE)
  }
  list(
    "unique" = p_unique, "item_type" = item_type,
    "strategies" = strategies
  )
}

match_pattern <- function(patterns, mat) {
  apply(patterns, 1, function(p) {
    idx <- which(apply(p == t(mat), 2, all))
    if (length(idx) == 0) idx <- NA
    idx
  })
}
