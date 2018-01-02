#' @return a matrix with two columns (Bayes factor and SE of approximation) and three rows:
#'     \itemize{
#'       \item \code{`bf_0u`}:  constrained vs. unconstrained (saturated) model
#'       \item \code{`bf_u0`}:  unconstrained vs. constrained model
#'       \item \code{`bf_00'`}: constrained vs. complement of inequality-constrained model
#'              (e.g., pi>.2 becomes pi<=.2; this assumes identical equality constraints for both models)
#'     }
#'
