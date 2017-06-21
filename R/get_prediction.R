
#' Get Predictions of Common Strategies
#'
#' Returns a vector of predictions for well-known choice strategies (TTB, WADD)
#' @param cueA cue values of Option A (-1/+1 = negative/positive; 0 = missing)
#' @param cueB cue values of Option B (see \code{cueA})
#' @param strategy strategy label, e.g., \code{"TTB"}, \code{"WADD"}, or \code{"WADDprob"}. See details.
#' @export
get_prediction <- function (cueA, cueB, strategy){
  pred <- switch(strategy,
                 "TTB" = {},
                 "EQW" = {},
                 "WADD" = {},
                 "GUESS" = {}
  )
  pred
}
