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
#'   \item{\code{"W1"}: }{Weak row/subject independence: Persons can be ordered on
#'                       an ordinal scale independent of items.}
#'   \item{\code{"W2"}: }{Weak column/item independence: Items can be ordered on
#'                       an ordinal scale independent of persons}
#'   \item{\code{"DC"}: }{Double cancellation: A necessary condition for a joint
#'                      ordering of (person,item) pairs and an additive representation
#'                      (i.e., an interval scale).}
#' }
#'
#' Note that axioms W1 and W2 jointly define the ISOP model by Scheiblechner
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
nirt_to_Ab <- function(N, M, options = 2, axioms = c("W1", "W2")){

  if (options == 2){
    idx <- outer(paste0("p", 1:N),
                 paste0("i", 1:M), paste, sep = ",")
    A <- matrix(0, 0, M*N, dimnames=list(NULL, idx))

    # weak row independence (participants ordered): (N-1)*M constraints
    if ("W1" %in% axioms){
      W1_tmp <- diag(1, N-1, N) - cbind(0, diag(N-1))
      for (m in 1:M){
        W1 <- cbind(matrix(0, N-1, (m-1)*N), W1_tmp, matrix(0, N-1, (M-m)*N))
        rownames(W1) <- rep("W1", nrow(W1))
        A <- rbind(A, W1 = W1)
      }
    }

    # weak column independence (participants ordered): N*(M-1) constraints
    if ("W2" %in% axioms){
      for (m in 1:(M-1)){
        for (n in 1:N){
          tmp <-  rep(0, ncol(A))
          tmp[(m-1)*N + c(n, n + N)] <- c(1, -1)
          A <- rbind(A, "W2" = tmp)
        }
      }
    }
    b <- rep(0, nrow(A))
    pt <- list("A" = A, "b" = b)

    if ("DC" %in% axioms){
      if (M < 3 | N < 3)
        stop("To test double cancellation (DC), at least a 3x3 table is required.")
      warning("Not yet implemented. Requires unions of convex polytopes.")

      # for a 3x3 table:
      # p12 < p21 & p23 < p32  =>  p13 < p31
      # <=> not(p12 < p21 & p23 < p32  & not(p13 < p31) )
      # <=> not(p12 < p21 & p23 < p32  & p13 > p31)
      #          4     2     8     6      7     3   [idx]
      DC.template <- matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0,
                              0, 0, 1, 0, 0, 0, -1, 0, 0,
                              0, 0, 0, 0, 0, -1, 0, 1, 0), 3, 9, byrow = TRUE)
      # similarly: p12 > p21 & p23 > p32  =>  p13 > p31
      #            => opposite: DC2 <- - DC1

      cols <- matrix(combn(M, 3), 3)
      rows <- matrix(combn(N, 3), 3)
      exclude <- vector("list", 2 * choose(M, 3)*choose(N,3))  # number of 3x3 tables
      cnt <- 1
      indices <- matrix(1:(N*M), N)   # indices in Ab representation
      for(j in 1:ncol(cols)){
        for(i in 1:ncol(rows)){
          A.tmp <- matrix(0, 3, M*N)
          A.idx <- indices[rows[,i], cols[,j]]
          A.tmp[,c(A.idx)] <- DC.template
          exclude[[cnt]]$A <- A.tmp
          exclude[[cnt + 1]]$A <- - A.tmp
          exclude[[cnt]]$b <- exclude[[cnt + 1]]$b <- rep(0, 3)
          cnt <- cnt + 2
        }
      }
      pt$exclude <- exclude
      # pt$b.exclude <- replicate(length(A.exclude), rep(0, 3), simplify = FALSE)
    }

  } else {
    # polytomous
    idx <-
      outer(paste0("p", 1:N),
            outer(paste0("i", 1:M),
                  paste0("o", 1:options - 1), paste, sep = ","), paste, sep = ",")
    coln <- c()
    for(i in 1:M) coln <- c(coln, t(idx[,i,-options]))
    A <- matrix(0, 0, M*N*(options-1), dimnames=list(NULL, coln))
    zeros <- matrix(0, options-1, ncol(A))
    ij_a <- ij_b <- matrix(0, options-1, options-1)
    ij_a[lower.tri(ij_a, TRUE)] <- -1
    ij_b[lower.tri(ij_b, TRUE)] <- 1   # P(t | ij_a) < P(t | ij_b)  for all categories

    # weak row independence (participants ordered): (N-1)*M*(options-1) constraints
    if ("W1" %in% axioms){
      for(n in 1:(N-1)){
        for (m in 1:M){
          idx_a <- (m-1)*N*(options-1)  + (n-1)*(options-1) + 1:(options-1)
          idx_b <- (m-1)*N*(options-1)  +     n*(options-1) + 1:(options-1)
          W1 <- zeros
          W1[,idx_a] <- ij_a
          W1[,idx_b] <- ij_b
          p_a <- strsplit(colnames(A)[idx_a], ",")[[1]][1]
          rownames(W1) <- paste0("W1_", p_a,",", colnames(A)[idx_b])
          A <- rbind(A, W1)
        }
      }
    }

    # weak column independence (items ordered): N*(M-1)*(options-1) constraints
    if ("W2" %in% axioms){
      for(n in 1:N){
        for (m in 1:(M-1)){
          idx_x <- (m-1)*N*(options-1)  + (n-1)*(options-1) + 1:(options-1)
          idx_y <-     m*N*(options-1)  + (n-1)*(options-1) + 1:(options-1)
          W2 <- zeros
          W2[,idx_x] <- ij_a
          W2[,idx_y] <- ij_b
          i_x <- strsplit(colnames(A)[idx_x], ",")[[1]][2]
          rownames(W2) <- paste0("W2_", i_x,",", colnames(A)[idx_y])
          A <- rbind(A, W2)
        }
      }
    }
    b <- rep(0, nrow(A))
    pt <- list("A" = A, "b" = b)
  }

  ################################ polytomous:
  # data: I x J x K contingency table
  # (unique persons, unique items, response categories)
  # possible: I <= N  , j<= M
  pt
}

