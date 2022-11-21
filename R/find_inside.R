#' Find a Point/Parameter Vector Within a Convex Polytope
#'
#' Finds the center/a random point that is within the convex polytope defined by the
#' linear inequalities \code{A*x <= b}  or by the convex hull over the vertices in the matrix \code{V}.
#'
#' @inheritParams count_binom
#' @param probs only for \code{A*x<b} representation: whether to add
#'     inequality constraints that the variables are probabilities (nonnegative and
#'     sum to 1 within each option)
#' @param options optional: number of options per item type (only for \eqn{A x \leq b} representation).
#'     Necessary to account for sum-to-one constraints within multinomial
#'     distributions (e.g., p_1 + p_2 + p_3 <= 1).
#'     By default, parameters are assumed to be independent.
#' @param random if \code{TRUE}, random starting values in the interior are generated.
#'     If \code{FALSE}, the center of the polytope is computed using linear programming.
#' @param boundary constant value \eqn{c} that is subtracted on the right-hand side
#'     of the order constraints, \eqn{A x \leq b - c}. This ensuresa that the
#'     resulting point is in the interior of the polytope and
#'     not at the boundary, which is important for MCMC sampling.
#'
#' @details
#' If vertices \code{V} are provided, a convex combination of the vertices is returned.
#' If \code{random=TRUE}, the weights are drawn uniformly from a Dirichlet distribution.
#'
#' If inequalities are provided via \code{A} and \code{b}, linear programming (LP) is used
#' to find the Chebyshev center of the polytope (i.e., the center of the largest ball that
#' fits into the polytope; the solution may not be unique). If \code{random=TRUE},
#' LP is used to find a random point (not uniformly sampled!) in the convex polytope.
#'
#' @examples
#' # inequality representation (A*x <= b)
#' A <- matrix(c(1, -1,  0, 1,  0,
#'               0,  0, -1, 0,  1,
#'               0,  0,  0, 1, -1,
#'               1,  1,  1, 1,  0,
#'               1,  1,  1, 0,  0,
#'               -1, 0, 0, 0, 0),
#'             ncol = 5, byrow = TRUE)
#' b <- c(0.5, 0, 0, .7, .4, -.2)
#' find_inside(A, b)
#' find_inside(A, b, random = TRUE)
#'
#'
#' # vertex representation
#' V <- matrix(c(
#'   # strict weak orders
#'   0, 1, 0, 1, 0, 1,  # a < b < c
#'   1, 0, 0, 1, 0, 1,  # b < a < c
#'   0, 1, 0, 1, 1, 0,  # a < c < b
#'   0, 1, 1, 0, 1, 0,  # c < a < b
#'   1, 0, 1, 0, 1, 0,  # c < b < a
#'   1, 0, 1, 0, 0, 1,  # b < c < a
#'
#'   0, 0, 0, 1, 0, 1,  # a ~ b < c
#'   0, 1, 0, 0, 1, 0,  # a ~ c < b
#'   1, 0, 1, 0, 0, 0,  # c ~ b < a
#'   0, 1, 0, 1, 0, 0,  # a < b ~ c
#'   1, 0, 0, 0, 0, 1,  # b < a ~ c
#'   0, 0, 1, 0, 1, 0,  # c < a ~ b
#'
#'   0, 0, 0, 0, 0, 0   # a ~ b ~ c
#' ), byrow = TRUE, ncol = 6)
#' find_inside(V = V)
#' find_inside(V = V, random = TRUE)
#' @export
find_inside <- function(A,
                        b,
                        V,
                        options = NULL,
                        random = FALSE,
                        probs = TRUE,
                        boundary = 1e-5){

  # (1) V-representation: convex combination of vertices
  if (!missing(V) && !is.null(V)){
    check_V(V)
    if (random){
      u <- c(rpdirichlet(1, rep(1, nrow(V)), nrow(V), drop_fixed = FALSE))
    } else {
      u <- rep(1/nrow(V), nrow(V))
    }
    p <- colSums(V * u)
    if (!inside(p, V = V))
      stop("No point found inside of convex hull of V.")


  # (2) Ab-representation: solution of quadratic program
  } else {
    if (missing(options) || is.null(options))
      options <- rep(2, ncol(A))
    if (probs){
      zeros <- rep(0, sum(options))
      check_Abokprior(A, b, options, zeros, zeros)
      tmp <- Ab_multinom(options, A, b, nonneg = TRUE)
      A <- tmp$A
      b <- tmp$b
    }

    if (random){
      u <- runif(ncol(A), 0, 1)
      B <- diag(ncol(A)) #matrix(runif(ncol(A)^2, -1,1), ncol(A))
      try (p <- quadprog::solve.QP(Dmat = B,
                                   dvec = u,
                                   Amat = - t(A),
                                   bvec = - b + boundary)$solution)
      p
    } else {
      # find analytic center of polytope:
      # https://math.stackexchange.com/questions/1377209/analytic-center-of-convex-polytope
      # http://stanford.edu/class/ee364a/lectures/problems.pdf  (page: 4-19)
      obj <- c(1, rep(0, ncol(A)))  # variables: (r, p1,p2,...)
      a <- apply(A, 1, norm, type = "2")
      Ar <- cbind(a, A)
      # QP much faster than LP, almost identical results:
      try (p <- solve.QP(Dmat = 1e-5*diag(ncol(Ar)),
                         dvec = obj,
                         Amat = - t(Ar),
                         bvec = - b + boundary)$solution[-1])
    }

    if (is.null(p) || !inside(p, A, b))
      stop ("No point found inside of A*x <=b. Maybe not all inequalities can be satisfied?!")
  }
  p
}


####### different LPs
# Rglpk
# ' @param tm_limit time limit for linear program in ms, see \link[Rglpk]{Rglpk_solve_LP}.
# dir <- rep("<=", nrow(A))
# p <- Rglpk::Rglpk_solve_LP(obj, Ar, dir, b, max = TRUE,
#                            control = list(tm_limit = 0))$solution[-1]
# lpSolve::lp("max", objective.in = obj, const.mat = Ar,
#             const.dir = "<=", const.rhs = b)$solution[-1]
