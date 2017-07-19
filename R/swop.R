
#' Strict Weak Order Polytope for 5 Elements and Ternary Choices
#'
#' Facet-defining inequalities of the strict weak order mixture model
#' for the 20 paired comparisons of 5 choice elements {a,b,c,d,e}
#' (Regenwetter & Davis-Stober, 2012).
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{\code{A}: }{Matrix with inequality constraints that define a polytope via A*x <= b.}
#'   \item{\code{b}: }{vector with upper bounds for the inequalities.}
#'   \item{\code{start}: }{A point in the polytope.}
#'   \item{\code{options}: }{A vector with the number of options (=3) per item type.}
#' }
#' @template ref_regenwetter2012
#' @examples
#' data(swop5)
#' tail(swop5$A)  # A*x <= b
#' tail(swop5$b)
#' swop5$start    # inside SWOP polytope
#' swop5$options  # number of options per item
#'
#' # check whether point is in polytope:
#' inside(swop5$start, swop5$A, swop5$b)
#'
#' \dontrun{
#' # get prior samples:
#' p <- sampling_multinomial(0, swop5$options,
#'                           swop5$A, swop5$b,
#'                           M = 1000, start = swop5$start)
#' colMeans(p)
#' apply(p[,1:5], 2, plot, type = "l")
#' }
"swop5"
