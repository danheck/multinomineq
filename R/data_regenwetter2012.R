#' Data: Ternary Risky Choices (Regenwetter & Davis-Stober, 2012)
#'
#' Raw data with choice frequencies for all 20 paired comparison of 5 gambles a, b, c, d, and e.
#' Participants could either choose "Option 1", "Option 2", or "indifferent" (ternary choice).
#' Each paired comparison (e.g., a vs. b) was repeated 45 times per participant.
#' The data include 3 different gamble sets and aimed at testing whether people
#' have transitive preferences (see Regenwetter & Davis-Stober, 2012).
#'
#' @format A matrix with 22 columns:
#' \describe{
#'   \item{\code{participant}: }{Participant number}
#'   \item{\code{gamble_set}: }{Gamble set}
#'   \item{\code{a>b}: }{Number of times a preferred over b}
#'   \item{\code{b>a}: }{Number of times b preferred over a}
#'   \item{\code{a=b}: }{Number of times being indifferent between a and b}
#' }
#' @template ref_regenwetter2012
#' @seealso The substantive model of interest was the strict weak order polytope (see \code{\link{swop5}}).
#'
#' @examples
#' data(regenwetter2012)
#' head(regenwetter2012)
#'
#' # check transitive preferences: strict weak order polytope (SWOP)
#' data(swop5)
#' tail(swop5$A, 3)
#' # participant 1, gamble set 1:
#' p1 <- regenwetter2012[1,-c(1:2)]
#' inside_multinom(p1, swop5$options, swop5$A, swop5$b)
#'
#' \dontrun{
#' # posterior samples
#' p <- sampling_multinom(regenwetter2012[1,-c(1:2)],
#'                        swop5$options, swop5$A, swop5$b,
#'                        M=1000, start = swop5$start)
#' colMeans(p)
#' apply(p[,1:6], 2, plot, type = "l")
#' ppp_multinom(p, p1, swop5$options)
#'
#' # Bayes factor
#' bf_multinom(regenwetter2012[1,-c(1:2)], swop5$options,
#'             swop5$A, swop5$b,
#'             M = 1000, cmin = 1, steps = seq(2000,75000,2000))
#' }
"regenwetter2012"

