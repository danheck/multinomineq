
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
#' }
#' @template ref_regenwetter2012
#' @examples
#' data(swop5)
#' tail(swop5$A)
#' tail(swop5$b)
#' swop5$start
"swop5"
