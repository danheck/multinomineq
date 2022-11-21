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
#'     (if provided, \code{A} and \code{b} are ignored).
#'     Similar as for \code{A}, columns of \code{V} omit the last value for each
#'     multinomial condition (e.g., a1,a2,a3,b1,b2 becomes a1,a2,b1).
#'     Note that this method is comparatively slow since it solves linear-programming problems
#'     to test whether a point is inside  a polytope (Fukuda, 2004) or to run the Gibbs sampler.
#'
#' @seealso \code{\link{Ab_to_V}} and \code{\link{V_to_Ab}} to change between A/b and V representation.
#' @examples
#' # linear order constraints:  x1<x2<x3<.5
#' A <- matrix(c(
#'   1, -1, 0,
#'   0, 1, -1,
#'   0, 0, 1
#' ), ncol = 3, byrow = TRUE)
#' b <- c(0, 0, .50)
#'
#' # vertices: admissible points (corners of polytope)
#' V <- matrix(c(
#'   0, 0, 0,
#'   0, 0, .5,
#'   0, .5, .5,
#'   .5, .5, .5
#' ), ncol = 3, byrow = TRUE)
#'
#' xin <- c(.1, .2, .45) # inside
#' inside(xin, A, b)
#' inside(xin, V = V)
#'
#' xout <- c(.4, .1, .55) # outside
#' inside(xout, A, b)
#' inside(xout, V = V)
#' @export
inside <- function(x, A, b, V) {
  if (!missing(V) && !is.null(V)) {
    check_Vx(V, x)
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
#' @inheritParams inside
#' @param k choice frequencies.
#'    For \code{inside_binom}: per item type (e.g.: a1,b1,c1,..)
#'    For \code{inside_multinom}: for all choice options ordered by item type
#'    (e.g., for ternary choices: a1,a2,a3, b1,b2,b3,..)
#' @param n only for \code{inside_binom}: number of choices per item type.
#' @param options only for \code{inside_multinom}: number of response options per item type.
#' @examples
#' ############ binomial
#' # x1<x2<x3<.50:
#' A <- matrix(c(
#'   1, -1, 0,
#'   0, 1, -1,
#'   0, 0, 1
#' ), ncol = 3, byrow = TRUE)
#' b <- c(0, 0, .50)
#' k <- c(0, 1, 5)
#' n <- c(10, 10, 10)
#' inside_binom(k, n, A, b)
#'
#' ############ multinomial
#' # two ternary choices:
#' #     (a1,a2,a3,   b1,b2,b3)
#' k <- c(1, 4, 10, 5, 9, 1)
#' options <- c(3, 3)
#' # a1<b1, a2<b2, no constraints on a3, b3
#' A <- matrix(c(
#'   1, -1, 0, 0,
#'   0, 0, 1, -1
#' ), ncol = 4, byrow = TRUE)
#' b <- c(0, 0)
#' inside_multinom(k, options, A, b)
#'
#' # V-representation:
#' V <- matrix(c(
#'   0, 0, 0, 0,
#'   0, 0, 0, 1,
#'   0, 1, 0, 0,
#'   0, 0, 1, 1,
#'   0, 1, 0, 1,
#'   1, 1, 0, 0,
#'   0, 1, 1, 1,
#'   1, 1, 0, 1,
#'   1, 1, 1, 1
#' ), 9, 4, byrow = TRUE)
#' inside_multinom(k, options, V = V)
#' @seealso \code{\link{inside}}
#' @export
inside_binom <- function(k, n, A, b, V) {
  check_kn(k, n)
  if (!is.null(dim(k)) && length(dim(k)) == 2) {
    x <- t(t(k) / n)
  } else {
    x <- k / n
  }
  if (missing(V) || is.null(V)) {
    inside(x, A, b)
  } else {
    inside(x, V = V)
  }
}

#' @rdname inside_binom
#' @export
inside_multinom <- function(k, options, A, b, V) {
  k_free <- drop_fixed(k, options)
  sel <- rep(1:length(options), options - 1)
  if (!is.null(dim(k_free)) && length(dim(k_free)) == 2) {
    n <- tapply(t(k), rep(1:length(options), options), sum)[sel]
    x <- t(t(k_free) / c(n))
  } else {
    n <- tapply(k, rep(1:length(options), options), sum)[sel]
    x <- k_free / n
  }
  if (missing(V) || is.null(V)) {
    inside(x, A, b)
  } else {
    inside(x, V = V)
  }
}

inside_V <- function(x, V, return_glpk = FALSE) {
  if (!is.null(dim(x)) && length(dim(x)) == 2) {
    return(apply(x, 1, inside_V, V = V, return_glpk = return_glpk))
  } else {
    # mb <- microbenchmark::microbenchmark(times = 20, fukuda = {

    # Fukuda 2004:  (similar to Smeulders 3.2)
    npar <- length(x) + 1
    obj <- c(-1, x) # z0, z
    mat <- cbind(-1, rbind(V, x)) # z0, z
    dir <- rep("<=", nrow(mat))
    rhs <- c(rep(0, nrow(V)), 1)
    bnd <- list(
      lower = list(ind = 1:npar, val = rep(-Inf, npar)),
      upper = list(ind = 1:npar, val = rep(Inf, npar))
    )
    glpk <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE, bounds = bnd)

    # }, smeulders = {
    #   # [CURRENTLY WRONG] Smeulders 2018: linear problem in 3.1
    # npar <- 1 + nrow(V)   # c(z, x_m)  where x_m is the mixture weight for row m of V
    # obj <- c(1, rep(0, npar - 1))  # minimize z
    # mat <- rbind(cbind(1, t(V)),            # sum x_m + z >= p
    #              c(0, rep(1, npar - 1)))   # sum -x_n >= -1
    # dir <- c(rep(">=", nrow(mat) - 1), "==")
    # rhs <- c(x, 1)
    # glpk2 <- Rglpk_solve_LP(obj, mat, dir, rhs, max = FALSE) # default: x_m, z >= 0
    #   # sum(glpk2$solution[-1])
    # },{
    #   # [CURRENTLY WRONG] Smeulders 2018: linear problem in 3.2
    # npar <- length(x) + 1  # c, y_ij
    # obj <- c(-1, x)       # -c, p_ij
    # mat <- rbind(cbind(-1, V),            # sum_ij  - c <= 0
    #              c(0, rep(1, npar - 1)))  # sum y_ij <= 1
    # dir <- c(rep("<=", nrow(mat) - 1),
    #          "<=")
    # rhs <- c(rep(0, nrow(V)), 1)
    # bnd <- list(lower = list(ind = 1:npar, val = rep(0, npar)),
    #             upper = list(ind = 1:npar, val = rep(Inf, npar)))
    # glpk3 <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE, bounds = bnd)
    # })
    # print(mb)
    # Ab <- V_to_Ab(V)
    # cat("Fukuda: ", glpk$optimum, " . Smeulders: ", glpk2$optimum, ". A*x<b: ", inside(x, Ab$A, Ab$b), "\n")

    if (return_glpk) {
      glpk$inside <- !glpk$optimum > 0
      glpk
    } else {
      return(!glpk$optimum > 0) # >0  --> outside (= separating hyperplane exists)
    }
  }
}
