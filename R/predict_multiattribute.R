
#' Get Predictions for Multiattribute Decisions
#'
#' Returns a vector of predictions for choice strategies (e.g., TTB, WADD)
#'
#' @param cueA cue values of Option A (-1/+1 = negative/positive; 0 = missing).
#'        If a matrix is provided, each row defines one item type.
#' @param cueB cue values of Option B (see \code{cueA}).
#' @param v cue validities: probabilities that cues lead to correct decision.
#'        Must be of the same length as the number of cues.
#' @param strategy strategy label, e.g., \code{"TTB"}, \code{"WADD"}, or \code{"WADDprob"}. See details.
#'
#' @template details_prediction
#' @return Predidcted choice pattern are encoded as negative = Option A, positive = Option B, and 0 = guessing.
#'         Different integers are used if error probabilities are assumed to be different.
#' @examples
#' ca <- c(1, -1, -1, 1)
#' cb <- c(-1, 1, -1, -1)
#' v <- c(.9, .8, .7, .6)
#' strategies <-
#' predict_multiattribute(ca, cb, v, "TTBprob")
#' predict_multiattribute(ca, cb, v, "TTB")
#' predict_multiattribute(ca, cb, v, "WADDprob")
#' predict_multiattribute(ca, cb, v, "WADD")
#' predict_multiattribute(ca, cb, v, "EQW")
#' predict_multiattribute(ca, cb, v, "GUESS")
#' @export
predict_multiattribute <- function (cueA, cueB, v, strategy){
  check_cues(cueA, cueB, v)

  if (length(strategy) != 1){
    # multiple strategies
    pred <- lapply(strategy, function(s)
                   predict_multiattribute(cueA, cueB, v, s))
    names(pred) <- strategy

  } else if (is.matrix(cueA) || is.data.frame(cueA) || strategy == "baseline"){
    # multiple cues
    if (!is.matrix(v))
      v <- matrix(v, nrow(cueA), length(v), byrow = TRUE)
    if (strategy == "baseline"){
      warning ("The upper boundary of the parameters must be set to c=1 manually.")
      pred <- 1:nrow(cueA)
    } else {
      pred <- mapply(predict_multiattribute,
                     cueA = as.list(data.frame(t(cueA))),
                     cueB = as.list(data.frame(t(cueB))),
                     v = as.list(data.frame(t(v))),
                     MoreArgs = list(strategy = strategy))
    }
    names(pred) <- rownames(cueA)

    if (strategy %in% c("WADDprob", "TTBprob"))
      attr(pred, "ordered") <- TRUE
    else
      attr(pred, "ordered") <- NULL
  } else {
    pred <- switch(strategy,
                   "TTBprob" = {
                     o <- order(v, decreasing = TRUE)
                     diff <- cueB[o] - cueA[o]
                     d2 <- diff[diff != 0]
                     cnt <- which.max(diff != 0)
                     ifelse(length(d2) == 0, 0, cnt * sign(d2[1]))
                   },
                   "TTB" = {
                     cnt <- predict_multiattribute(cueA, cueB, v, "TTBprob")
                     sign(cnt)
                   },
                   "EQW" = {
                     cnt_diff <- sum(cueB == 1) - sum(cueA == 1)
                     sign(cnt_diff)
                   },
                   "WADDprob" = {
                     logodds <- log(v/(1-v))
                     sum(logodds[cueA == -1 & cueB == 1]) -
                       sum(logodds[cueA == 1 & cueB == -1])
                   },
                   "WADD" = {
                     ev <- predict_multiattribute(cueA, cueB, v, "WADDprob")
                     sign(ev)
                   },
                   "GUESS" = {
                     0.0
                   })
  }
  pred
}



#' Transform Vector of Predictions to Polytope
#'
#' Transforms ordered item-type predictions to polytope definition. This allows to use Monte-Carlo methods to compute the Bayes factor if the number of item types is large (\code{\link{compute_bf}}).
#'
#' @param prediction a numeric vector with probabilistic choice predictions.
#'   See details and \code{\link{predict_multiattribute}}.
#' @param c upper boundary of error probabilities
#'
#' @template details_prediction
#' @details
#' Computes the matrix A and the vector b that define a polytope via {x: A*x < c}.
#' Only works for models without guessing predictions and equality constraints (i.e., different error probabilities per item type)
#' @examples
#' k <- c(1,7,2,5)
#' n <- rep(10, 4)
#' # predict: A,A,B,B
#' pred <- c(-1,2,-3,4)
#' poly <- as_polytope(pred, c= .5)
#' cbind(poly$A, b = poly$b)
#'
#' bf <- exp(compute_marginal(k, n, pred) - compute_marginal(k, n, 1:4, c = 1))
#' bf
#' compute_bf(k, n, poly$A, poly$b, prior = c(1,1), M = 2e5)
#' @export
as_polytope <- function(prediction, c = .50){
  pu <- get_par_unique(prediction)
  npar <- get_par_number(prediction)
  if (any(prediction == 0) || length(unique(prediction)) != npar || npar < 2)
    stop ("Not working for guesses or if there are less parameters than item types.")

  # prediction negative: 0<P(b)<c ; positive: c<P(b)<1
  s <- sign(prediction)
  A <- rbind(diag(s), diag(-s))
  colnames(A) <- paste0("par", pu)
  b <- c(ifelse(s == -1, 0, 1),
         -s * c)

  # linear order constraints:
  Aloc <- matrix(0, npar - 1, npar)
  bloc <- rep(0, npar - 1)
  for (i in 2:npar){
    # find correct indices to add order constraints
    if (prediction[i] < 0){
      if (prediction[i - 1] < 0){
        # -/-   =   A/A
        if (abs(prediction[i - 1]) < abs(prediction[i]) )  # b1 < b2
          Aloc[i-1,c(i-1, i)] <- c(1, -1)
        else                                               # b1 > b2
          Aloc[i-1,c(i-1, i)] <- c(-1, 1)
        # bloc[i - 1] <- 0

      } else {
        # -/+    =  A/B
        if (abs(prediction[i - 1]) < abs(prediction[i]) ){  # b1  < 1-b2
          Aloc[i-1,c(i-1, i)] <- c(-1,-1)
          bloc[i - 1] <- -1
        }else{
          Aloc[i-1,c(i-1, i)] <- c(1,1)                    # b1  > 1-b2
          bloc[i - 1] <- 1
        }
      }

    } else {
      if (prediction[i - 1] > 0){
        # +/+   =   B/B
        if (abs(prediction[i - 1]) < abs(prediction[i]) )  # 1-b1 < 1-b2
          Aloc[i-1,c(i-1, i)] <- c(-1, 1)
        else
          Aloc[i-1,c(i-1, i)] <- c(1, -1)                  # 1-b1 > 1-b2
        # bloc[i - 1] <- 0

      } else {
        # +/-   = B/A
        if (abs(prediction[i - 1]) < abs(prediction[i]) ){   # 1-b1  < b2
          Aloc[i-1,c(i-1, i)] <- c(1,1)
          bloc[i - 1] <- 1
        }else{
          Aloc[i-1,c(i-1, i)] <- c(-1,-1)                    # 1-b1  > b2
          bloc[i - 1] <- -1
        }
      }
    }
  }
  list(A = rbind(A, Aloc), b = c(b, bloc))
}

# proB   <- estimate_par(k, n, - prediction, prob = FALSE)
# proA <- estimate_par(k, n,   prediction, prob = FALSE)
# k <- ifelse(sign(prediction) == -1, proB, proA)
# n_poly <- proA + proB
# k <- k[order(abs(prediction))]
# if (npar > 1 && c != 1){
#   Aloc <- diag(1, npar - 1, npar) - cbind(0, diag(npar - 1))
#   A <- rbind(A, Aloc)
#   b <- c(b, rep(0, npar - 1))
# }
# s1 <- stratsel:::encompassing_stepwise(k, n, poly$A, poly$b, c(1,1), 5e4, 9)
# s2 <- stratsel:::encompassing_stepwise(rep(0,4),rep(0,4), poly$A, poly$b, c(1,1), 5e4, 7)
# s1$int  /s2$int
# log(s1$int  /s2$int)
