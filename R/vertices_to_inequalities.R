#' Transform Vertex/Inequality Representations of Polytopes
#'
#' For convex polytopes: Uses \code{rPorta} (\url{https://github.com/TasCL/rPorta})
#' to transform the vertex representation to/from the inequality representation.
#'
#' @param V a matrix with one vertex of a polytope per row
#'        (i.e., the admissible preference orders of a random utility model)
#' @details
#' When defining a polytope via \code{A} and \code{b}, constraints are automatically
#' added to ensure that all parameters are probabilities:  0 <= x <= 1.
#' Note that the transformation can be very slow and might require days of computing.
#' @examples
#' \dontrun{
#' V <- matrix(c(0, 0, 0,
#'               1, 0, 0,
#'               0, .5, 0,
#'               0, 0, .5,
#'               .5, .5, .5,
#'               0, .5, 0), ncol = 3, byrow = TRUE)
#' vertex_to_ineq(V)
#'
#' A <- matrix(c(1, -1, 0,   # x1 < x2
#'               0, 1, -1),  # x2 < x3
#'             ncol = 3, byrow = TRUE)
#' b <- c(0, 0)
#' ineq_to_vertex(A, b)
#' }
#' @export
vertex_to_ineq <- function (V){
  if (!requireNamespace("rPorta", quietly = TRUE))
    stop ("The pacakge 'rPorta' is required (https://github.com/TasCL/rPorta).",
            call. = FALSE)
  check_V(V)
  # use rPorta
  poi <- rPorta::as.poiFile(V)
  ieq <- rPorta::traf(poi)
  unlink("porta.log")
  if (!all(ieq@inequalities@sign == -1))
    warning ("Inequalities are not '<=' (i.e., Porta::ieq sign != -1)")
  Ab <- ieq@inequalities@num / ieq@inequalities@den
  list("A" = Ab[,1:(ncol(Ab)-1 )], "b" = Ab[,ncol(Ab)])
}


#' @inheritParams compute_bf
#' @rdname vertex_to_ineq
#' @export
ineq_to_vertex <- function (A, b){
  if (!requireNamespace("rPorta", quietly = TRUE))
    stop ("The pacakge 'rPorta' is required (https://github.com/TasCL/rPorta).",
          call. = FALSE)
  check_Ab(A, b)
  S <- ncol(A)
  A <- rbind(A, diag(S), -diag(S))
  b <- c(b, rep(1, S), rep(0, S))

  ieq <- rPorta::as.ieqFile(cbind(A, b), sign = rep(- 1, length(b)))
  poi <- rPorta::traf(ieq)
  unlink("porta.log")
  V <- poi@convex_hull@num / poi@convex_hull@den
  V
}
