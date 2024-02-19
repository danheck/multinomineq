#' @return
#' a \code{strategy} object (a list) with the entries:
#' \describe{
#' \item{\code{pattern}: }{a numeric vector encoding the predicted choice pattern by the sign
#'        (negative = Option A, positive = Option B, 0 = guessing).
#'       Identical error probabilities are encoded by using the same absolute number
#'       (e.g., \code{c(-1,1,1)} defines one error probability with A,B,B predictions).}
#' \item{\code{c}: }{upper boundary of error probabilities}
#' \item{\code{ordered}: }{whether error probabilities are linearly ordered by their absolute value in \code{pattern} (largest error: smallest absolute number)}
#' \item{\code{prior}: }{a numeric vector with two positive values specifying the shape parameters of the beta prior distribution (truncated to the interval \code{[0,c]}}
#' \item{\code{label}: }{strategy label (optional)}
#' }
#'
