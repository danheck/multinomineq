#' Nonparametric Item Response Theory (NIRT)
#'
#' Provides the inequality constraints on choice probabilities implied by
#' nonparametric item response theory (NIRT; Karabatsos, 2001).
#'
#' @param N number of persons / rows in item-response table
#' @param M number of items / columns in item-response table
#' @param options number of item categories/response options. If \code{options=2},
#'    a dichotomous NIRT for product-binomial data is returned.
#' @param axioms which axioms should be included in the polytope representation \eqn{A*x <= b}?
#'     See details.
#'
#' @details
#' In contrast to parametric IRT models (e.g., the 1-parameter-logistic Rasch model),
#' NIRT does not assume specific parametric shapes of the item-response and person-response
#' functions. Instead, the necessary axioms for a unidimensional representation of
#' the latent trait are tested directly.
#'
#' The axioms are as follows:
#' \itemize{
#'   \item{\code{"WI1"}}{Weak row/subject independence: Persons can be ordered on
#'                       an ordinal scale independent of items.}
#'   \item{\code{"WI2"}}{Weak column/item independence: Items can be ordered on
#'                       an ordinal scale independent of persons}
#'   \item{\code{"DC"}}{Double cancellation: A necessary condition for a joint
#'                      ordering of (person,item) pairs and an additive representation
#'                      (i.e., an interval scale).}
#' }
#'
#' Note that axioms WI1 and WI2 jointly define the ISOP model by Scheiblechner
#' (1995; isotonic ordinal probabilistic model) and the double homogeneity model
#' by Mokken (1971). If DC is added, we obtain the ADISOP model
#' by Scheiblechner (1999; ).
#'
#' @examples
#' # 5 persons, 3 items
#' nirt_to_Ab(5, 3)
#' @template ref_karabatsos2001
#' @template ref_karabatsos2004
#' @template ref_mokken1971
#' @template ref_scheiblechner1995
#' @template ref_scheiblechner1999
#' @export
nirt_to_Ab <- function(N, M, options = 2, axioms = c("WI1", "WI2")){

  idx <- outer(paste0("p", 1:N),
               paste0("i", 1:M), paste, sep = ",")
  A <- matrix(0, 0, M*N, dimnames=list(NULL, idx))

  # weak row independence (participants ordered): (N-1)*M constraints
  if ("WI1" %in% axioms){
    WI1_tmp <- diag(1, N-1, N) - cbind(0, diag(N-1))
    for (m in 1:M){
      WI1 <- cbind(matrix(0, N-1, (m-1)*N), WI1_tmp, matrix(0, N-1, (M-m)*N))
      rownames(WI1) <- rep("WI1", nrow(WI1))
      A <- rbind(A, WI1 = WI1)
    }
  }

  # weak column independence (participants ordered): N*(M-1) constraints
  if ("WI2" %in% axioms){
    for (m in 1:(M-1)){
      for (n in 1:N){
        tmp <-  rep(0, ncol(A))
        tmp[(m-1)*N + c(n, n + N)] <- c(1, -1)
        A <- rbind(A, "WI2" = tmp)
      }
    }
  }

  ################################ polytomous:
  # data: I x J x K contingency table
  # (unique persons, unique items, response categories)
  # possible: I <= N  , j<= M


  b <- rep(0, nrow(A))
  list("A" = A, "b" = b)
}

