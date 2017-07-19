#' Check Whether Points are Inside a Convex Polytope
#'
#' Determines whether a point \code{x} is inside a convex poltyope by checking whether
#' (1) all inequalities \code{A*x <= b} are satisfied or
#' (2) the point \code{x} is in the convex hull of the vertices in \code{V}.
#'
#' @param x a vector of length equal to the number of columns of \code{A} or \code{V}
#' (i.e., a single point in D-dimensional space) or matrix of points/vertices (one per row).
#' @param A a matrix with one row for each linear inequality constraint and one
#'   column for each of the free parameters. The parameter space is defined
#'   as all probabilities \code{x} that fulfill the order constraints  \code{A*x <= b}.
#' @param b a vector of the same length as the number of rows of \code{A}.
#' @param V a matrix of vertices (one per row) that define the polytope of
#'     admissible parameters as the convex hull over these points
#'     (if provided, \code{A} and \code{b} are ignored). Note that this method
#'     is comparatively slow since it solves a linear-programming problem for each point (Fukuda, 2004).
#'
#' @seealso \code{\link{ineq_to_vertex}} to change between A/b and V representation.
#' @examples
#' # linear order constraints:  x1<x2<x3<.5
#' A <- matrix(c(1,-1, 0,
#'               0, 1,-1,
#'               0, 0, 1), ncol = 3, byrow = TRUE)
#' b <- c(0, 0, .50)
#' # vertices: admissible points (corners of polytope)
#' V <- matrix(c( 0, 0, 0,
#'                0, 0,.5,
#'                0,.5,.5,
#'               .5,.5,.5), ncol = 3, byrow = TRUE)
#'
#' xin <- c(.1, .2, .45)  # inside
#' inside(xin, A, b)
#' inside(xin, V = V)
#'
#' xout <- c(.4, .1, .55)  # outside
#' inside(xout, A, b)
#' inside(xout, V = V)
#' @export
inside <- function(x, A, b, V){
  if (!missing(V) && !is.null(V)){
    i <- inside_V(x, V)
  } else {
    check_Abx(A, b, x)
    i <- inside_Ab(matrix(x, ncol = ncol(A)), A, b)
  }
  as.logical(c(i))
}

#' Check Whether Choice Frequencies are in Polytope
#'
#' Computes relative choice frequencies and checks whether these are in the polytope defined
#' via (1) \code{A*x <= b} or (2) by the convex hull of a set of vertices \code{V}.
#'
#' @param k choice frequencies.
#'    For \code{inside_binomial}: per item type (e.g.: a1,b1,c1,..)
#'    For \code{inside_multinomial}: for all choice options ordered by item type
#'    (e.g., for ternary choices: a1,a2,a3, b1,b2,b3,..)
#' @param n only for \code{inside_binomial}: number of choices per item type.
#' @param options only for \code{inside_multinomial}: number of response options per item type.
#' @examples
#' ############ binomial
#' # x1<x2<x3<.50:
#' A <- matrix(c(1,-1,0,
#'               0,1,-1,
#'               0,0, 1), ncol=3, byrow=TRUE)
#' b <- c(0, 0, .50)
#' k <- c( 0, 1, 5)
#' n <- c(10,10,10)
#' inside_binomial(k, n, A, b)
#'
#' ############ multinomial
#' # two ternary choices:
#' #     (a1,a2,a3,   b1,b2,b3)
#' k <- c(1,4,10,     5,9,1)
#' options <- c(3, 3)
#' # a1<b1, a2<b2:
#' A <- matrix(c(1,-1,0, 0,
#'               0, 0,1,-1), ncol=4, byrow=TRUE)
#' b <- c(0, 0)
#' inside_multinomial(k, options, A, b)
#' @seealso \code{\link{inside}}
#' @export
inside_binomial <- function(k, n, A, b, V){
  x <- k / n
  if (missing(V)){
    inside(x, A, b)
  } else {
    inside(x, V = V)
  }
}

#' @rdname inside_binomial
#' @export
inside_multinomial <- function(k, options, A, b, V){
  k_free <- k[- cumsum(options)]
  sel <- rep(1:length(options), options - 1)
  n <- tapply(k, rep(1:length(options), options), sum)[sel]
  x <- k_free / n
  if (missing(V)){
    inside(x, A, b)
  } else {
    inside(x, V = V)
  }
}

inside_V <- function (x, V){
  check_Vx(V, x)
  if (!is.null(dim(x))){
    return (apply(x, 1, inside_V, V = V))
  } else {
    npar <- length(x) + 1
    obj <- c(-1, x)  # z0, z
    mat <- cbind(-1, rbind(V, x))  # z0, z
    dir <- rep("<=", nrow(mat))
    rhs <- c(rep(0, nrow(V)), 1)
    bnd <- list(lower = list(ind = 1:npar, val = rep(-Inf, npar)),
                upper = list(ind = 1:npar, val = rep(Inf, npar)))
    glpk <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE, bounds = bnd)
    return (!glpk$optimum > 0)  # >0  --> outside (= separating hyperplane exists)
  }
}
