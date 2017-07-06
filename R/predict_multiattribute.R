
#' Strategy Predictions for Multiattribute Decisions
#'
#' Returns a list defining the predictions of different choice strategies (e.g., TTB, WADD)
#'
#' @param cueA cue values of Option A (-1/+1 = negative/positive; 0 = missing).
#'        If a matrix is provided, each row defines one item type.
#' @param cueB cue values of Option B (see \code{cueA}).
#' @param v cue validities: probabilities that cues lead to correct decision.
#'        Must be of the same length as the number of cues.
#' @param strategy strategy label, e.g., \code{"TTB"}, \code{"WADD"}, or \code{"WADDprob"}.
#'        Can be a vector. See details.
#'
#' @return
#' a \code{strategy} object (a list) with the entries:
#' \itemize{
#' \item{\code{pattern}: }{a numeric vector encoding the predicted choice pattern by the sign
#'        (negative = Option A, positive = Option B, 0 = guessing).
#'       Identical error probabilities are encoded by using the same absolute number
#'       (e.g., \code{c(-1,1,1)} defines one error probability with A,B,B predictions).}
#' \item{\code{c}}{upper boundary of error probabilities}
#' \item{\code{ordered}: }{whether error probabilities are linearly ordered by their absolute value in \code{pattern} (largest error: smallest absolute number)}
#' \item{\code{prior}: }{a numeric vector with two positive values specifying the shape parameters of the beta prior distribution (truncated to the interval \code{[0,c]}}
#' \item{\code{label}: }{strategy label}
#' }
#'
#' @examples
#' ca <- c(1, -1, -1, 1)
#' cb <- c(-1, 1, -1, -1)
#' v <- c(.9, .8, .7, .6)
#'
#' predict_multiattribute(ca, cb, v, "TTB")
#' predict_multiattribute(ca, cb, v, "WADDprob")
#'
#' predict_multiattribute(rbind(c1 = ca, c2 = c(1,1,1,1)),
#'                        rbind(cb, c(1,1,-1,1)),
#'                        v, "WADDprob")
#' @export
predict_multiattribute <- function (cueA, cueB, v, strategy,
                                    c= .50, prior = c(1, 1)){
  check_cues(cueA, cueB, v)

  if (length(strategy) != 1){
    # multiple strategies
    strat.list <- lapply(strategy, function(s)
      predict_multiattribute(cueA, cueB, v, s, c, prior))
    names(strat.list) <- strategy
    return(strat.list)

  } else if (!is.null(dim(cueA)) || strategy == "baseline"){
    # multiple cues
    if (!is.matrix(v))
      v <- matrix(v, nrow(cueA), length(v), byrow = TRUE)
    if (strategy == "baseline"){
      item.list <- list(pattern = 1:nrow(cueA),
                        c = 1, ordered = FALSE,
                        prior = prior,
                        label = strategy)
    } else {
      if (length(c) == 1) c <- rep(c, nrow(cueA))
      tmp <- mapply(predict_multiattribute,
                    cueA = as.list(data.frame(t(cueA))),
                    cueB = as.list(data.frame(t(cueB))),
                    v = as.list(data.frame(t(v))),
                    c = c,
                    MoreArgs = list(strategy = strategy, prior = prior),
                    SIMPLIFY = FALSE)
      item.list <- tmp[[1]]
      item.list$pattern <- sapply(tmp, "[[", "pattern")
    }
    names(item.list$pattern) <- rownames(cueA)
    return(item.list)

  } else {
    # single strategy
    pred <- switch(strategy,
                   "TTBprob" = {
                     o <- order(v, decreasing = TRUE)
                     diff <- cueB[o] - cueA[o]
                     d2 <- diff[diff != 0]
                     cnt <- which.max(diff != 0)
                     ifelse(length(d2) == 0, 0, cnt * sign(d2[1]))
                   },
                   "TTB" = {
                     cnt <- predict_multiattribute(cueA, cueB, v, "TTBprob")$pattern
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
                     ev <- predict_multiattribute(cueA, cueB, v, "WADDprob")$pattern
                     sign(ev)
                   },
                   "GUESS" = {
                     0.0
                   },
                   {stop("Strategy not supported.")})
    pred.list <- list(pattern = pred,
                      c = c,
                      ordered = grepl("prob", strategy),
                      prior = prior,
                      label = strategy)
    class(pred.list) <- "strategy"
    return(pred.list)
  }
}



#' Transform Pattern of Predictions to Polytope
#'
#' Transforms ordered item-type predictions to polytope definition.
#' This allows to use Monte-Carlo methods to compute the Bayes factor
#' if the number of item types is large (\code{\link{compute_bf}}).
#'
#' @param strategy a decision strategy returned by \code{\link{predict_multiattribute}}.
#' @param c upper boundary of error probabilities
#'
#' @details
#' Note: Only works for models without guessing predictions and
#'       without equality constraints (i.e., requires separate error probabilities per item type)
#'
#' @return a list containing the matrix \code{A} and the vector \code{b}
#'      that define a polytope via \code{A*x <= b}.
#' @examples
#' # define a strategy:
#' strat <- list(prediction = c(1,-2,3,4),  # B,A,B,B
#'               c = .5, ordered = TRUE,
#'               prior = c(1,1), label = "test")
#' as_polytope(strat)
#' @export
as_polytope <- function(strategy){
  check_strategy(strategy)
  pattern <- strategy$pattern
  pu <- get_error_unique(pattern)
  n_error <- get_error_number(pattern)
  if (any(pattern == 0) ||
      length(unique(pattern)) != n_error || n_error < 2)
    stop ("Not working for guesses or if there are less errors than item types.")

  # pattern negative: 0<P(b)<c ; positive: c<P(b)<1
  s <- sign(pattern)
  A <- rbind(diag(s), diag(-s))
  colnames(A) <- paste0("error", pu)
  b <- c(ifelse(s == -1, 0, 1), -s * c)

  if (ordered){
    # linear order constraints:
    Aloc <- matrix(0, n_error - 1, n_error)
    bloc <- rep(0, n_error - 1)
    for (i in 2:n_error){
      # find correct indices to add order constraints
      if (pattern[i] < 0){
        if (pattern[i - 1] < 0){
          # -/-   =   A/A
          if (abs(pattern[i - 1]) < abs(pattern[i]) )  # b1 < b2
            Aloc[i-1,c(i-1, i)] <- c(1, -1)
          else                                               # b1 > b2
            Aloc[i-1,c(i-1, i)] <- c(-1, 1)
          # bloc[i - 1] <- 0

        } else {
          # -/+    =  A/B
          if (abs(pattern[i - 1]) < abs(pattern[i]) ){  # b1  < 1-b2
            Aloc[i-1,c(i-1, i)] <- c(-1,-1)
            bloc[i - 1] <- -1
          }else{
            Aloc[i-1,c(i-1, i)] <- c(1,1)                    # b1  > 1-b2
            bloc[i - 1] <- 1
          }
        }

      } else {
        if (pattern[i - 1] > 0){
          # +/+   =   B/B
          if (abs(pattern[i - 1]) < abs(pattern[i]) )  # 1-b1 < 1-b2
            Aloc[i-1,c(i-1, i)] <- c(-1, 1)
          else
            Aloc[i-1,c(i-1, i)] <- c(1, -1)                  # 1-b1 > 1-b2
          # bloc[i - 1] <- 0

        } else {
          # +/-   = B/A
          if (abs(pattern[i - 1]) < abs(pattern[i]) ){   # 1-b1  < b2
            Aloc[i-1,c(i-1, i)] <- c(1,1)
            bloc[i - 1] <- 1
          }else{
            Aloc[i-1,c(i-1, i)] <- c(-1,-1)                    # 1-b1  > b2
            bloc[i - 1] <- -1
          }
        }
      }
    }
    A <- rbind(A, Aloc)
    b <- c(b, bloc)
  }
  list(A = A, b = b)
}

# proB   <- estimate_error(k, n, - pattern, prob = FALSE)
# proA <- estimate_error(k, n,   pattern, prob = FALSE)
# k <- ifelse(sign(pattern) == -1, proB, proA)
# n_poly <- proA + proB
# k <- k[order(abs(pattern))]
# if (n_error > 1 && c != 1){
#   Aloc <- diag(1, n_error - 1, n_error) - cbind(0, diag(n_error - 1))
#   A <- rbind(A, Aloc)
#   b <- c(b, rep(0, n_error - 1))
# }
# s1 <- stratsel:::encompassing_stepwise(k, n, poly$A, poly$b, c(1,1), 5e4, 9)
# s2 <- stratsel:::encompassing_stepwise(rep(0,4),rep(0,4), poly$A, poly$b, c(1,1), 5e4, 7)
# s1$int  /s2$int
# log(s1$int  /s2$int)
