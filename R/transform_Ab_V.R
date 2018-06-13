#' Transform Vertex/Inequality Representation of Polytope
#'
#' For convex polytopes: Requires \code{rPorta} (\url{https://github.com/TasCL/rPorta})
#' to transform the vertex representation to/from the inequality representation.
#'
#' @param V a matrix with one vertex of a polytope per row
#'        (e.g., the admissible preference orders of a random utility model or any other theory).
#'        Since the values have to sum up to one within each multinomial condition,
#'        the last value of each multinomial is omitted
#'        (e.g., the prediction 1-0-0/0-1 for a tri and binomial becomes 1-0/0).
#' @inheritParams count_multinom
#' @details
#' Choice models can be represented as polytopes if they assume a latent
#' mixture over a finite number preference patterns (random preference model).
#' For the general approach and theory underlying binary and ternary choice models,
#' see Regenwetter et al. (2012, 2014, 2017).
#'
#' Note that the transformation can be very slow and might require days or months
#' of computing or not converge at all!
#'
#' @template ref_regenwetter2012
#' @template ref_regenwetter2014
#' @template ref_regenwetter2017
#' @examples
#' \dontrun{
#' ######## (requires rPorta) ########
#'
#' ### binary choice:
#' # linear order: x1 < x2 < x3 < .50
#' # (cf. WADDprob in ?predict_multiattribute)
#' A <- matrix(c(1, -1,  0,
#'               0,  1, -1,
#'               0,  0,  1),
#'             ncol = 3, byrow = TRUE)
#' b <- c(0, 0, .50)
#' Ab_to_V(A, b)
#'
#'
#' ### binary choice polytope:
#' # choice options: {prefer_a, prefer_b}
#' # column order of vertices: (ab, ac, bc)
#' # with:  ij = 1  <=>  utility(i) > utility(j)
#' V <- matrix(c(1, 1, 1,  # c < b < a
#'               1, 1, 0,  # b < c < a
#'               0, 1, 1,  # c < a < b
#'               0, 0, 1   # a < c < b
#'             ), ncol = 3, byrow = TRUE)
#' V_to_Ab(V)
#'
#'
#' ### ternary choice (Regenwetter & Davis-Stober, 2012)
#' # choice options:  {prefer_a, indifferent, prefer_b}
#' # column order:    (ab,ba,  ac,ca,  bc,cb)
#' # with:            ij = 1  <=> utility(i) > utility(j)
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
#' V_to_Ab(V)
#' }
#' @export
V_to_Ab <- function (V){
  # options <- check_V(V, options)
  check_V(V)
  if (!requireNamespace("rPorta", quietly = TRUE))
    stop ("The pacakge 'rPorta' is required (https://github.com/TasCL/rPorta).",
          call. = FALSE)
  # use rPorta
  poi <- rPorta::as.poiFile(V)
  ieq <- rPorta::traf(poi)
  unlink("porta.log")
  if (!all(ieq@inequalities@sign == -1)){
    warning ("Inequalities are not '<=' (i.e., Porta::ieq sign != -1).",
             "\n  Complete Porta object is returned.")
    return (ieq)
  }
  Ab <- ieq@inequalities@num / ieq@inequalities@den
  A <- Ab[,1:(ncol(Ab)-1 )]
  colnames(A) <- colnames(V)
  list("A" = A, "b" = Ab[,ncol(Ab)])
}


#' @inheritParams count_multinom
#' @rdname V_to_Ab
#' @param options number of choice options per item type.
#'    Can be a vector \code{options=c(2,3,4)} if item types have 2/3/4 choice options.
#' @details
#' For binary choices (\code{options=2}), additional constraints are added to \code{A} and \code{b}
#' to ensure that all dimensions of the polytope satisfy:  0 <= p_i <= 1.
#' For ternary choices (\code{options=3}), constraints are added to ensure that 0 <= p_1+p_2 <=1
#' for pairwise columns (1+2, 3+4, 5+6, ...). See \code{\link{Ab_multinom}}.
#'
#' @export
Ab_to_V <- function (A, b, options = 2){
  if (!requireNamespace("rPorta", quietly = TRUE))
    stop ("The pacakge 'rPorta' is required (https://github.com/TasCL/rPorta).",
          call. = FALSE)

  options <- check_Ab(A, b, options)
  tmp <- Ab_multinom(options, A, b, nonneg = TRUE)
  A <- tmp$A
  b <- tmp$b
  ieq <- rPorta::as.ieqFile(cbind(A, b), sign = rep(- 1, length(b)))
  poi <- rPorta::traf(ieq)
  unlink("porta.log")
  V <- poi@convex_hull@num / poi@convex_hull@den
  colnames(V) <- colnames(A)
  V
}


#' Get Constraints for Product-Multinomial Probabilities
#'
#' Get or add inequality constraints (or vertices) to ensure that multinomial probabilities are
#' positive and sum to one for all choice options within each item type.
#'
#' @inheritParams count_multinom
#' @inheritParams Ab_to_V
#' @param nonneg whether to add constraints that probabilities must be nonnegative
#' @details
#' If \code{A} and \code{b} are provided, the constraints are added to these inequality constraints.
#' @seealso \code{\link{add_fixed}}
#'
#' @examples
#' # three binary and two ternary choices:
#' options <- c(2,2,2, 3,3)
#' Ab_multinom(options)
#' Ab_multinom(options, nonneg = TRUE)
#' @export
Ab_multinom <- function (options, A = NULL, b = NULL, nonneg = FALSE){
  S <- sum(options - 1)
  sum_to_one <- matrix(0, length(options), S)
  cnt <- 0
  for (i in 1:length(options)){
    idx <- seq.int(1, options[i] - 1) + cnt
    sum_to_one[i,idx] <- 1
    cnt <- cnt + options[i] - 1
  }
  if (nonneg){
    A_new <- rbind(A, - diag(S), sum_to_one)
    b_new <- c(b, rep(0, S), rep(1, length(options)))
  } else {
    A_new <- rbind(A, sum_to_one)
    b_new <- c(b, rep(1, length(options)))
  }
  list("A" = A_new, "b" = b_new)
}

# ' # convex hull of vertices (binomial and trinomial)
# ' V <- matrix(c(1,  0,0,
# '               0,  1,0,
# '               0, .5,.5), 3, 3, byrow = TRUE)
# ' V_multinom(options = c(2, 3),  V)
# ' @rdname Ab_multinom
# ' @inheritParams V_to_Ab
# ' @export
# V_multinom <- function (options, V){
#   see: add_fixed
# }
