# make list with default strategy options
as_strategy <- function(pattern, c = .50, ordered = TRUE, prior = c(1, 1)) {
  strategy <- list(pattern = pattern, c = c, ordered = ordered, prior = prior)
  class(strategy) <- "strategy"
  strategy
}


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
#' @param c defines the upper boundary for the error probabilities
#' @param prior defines the prior distribution for the error probabilities
#'          (i.e., truncated independent beta distributions \code{dbeta(prior[1], prior[2])} )
#'
#' @return
#' a \code{strategy} object (a list) with the entries:
#' \describe{
#' \item{\code{pattern}: }{a numeric vector encoding the predicted choice pattern by the sign
#'        (negative = Option A, positive = Option B, 0 = guessing).
#'       Identical error probabilities are encoded by using the same absolute number
#'       (e.g., \code{c(-1,1,1)} defines one error probability with A,B,B predictions).}
#' \item{\code{c}: }{upper boundary of error probabilities}
#' \item{\code{ordered}: }{whether error probabilities are linearly ordered by their absolute value in \code{pattern} (largest error: smallest absolute number)}
#' \item{\code{prior}: }{a numeric vector with two positive values specifying the shape parameters of the beta prior distribution (truncated to the interval \code{[0,c]}}
#' \item{\code{label}: }{strategy label}
#' }
#'
#' @examples
#' # single item type
#' v <- c(.9, .8, .7, .6)
#' ca <- c(1, -1, -1, 1)
#' cb <- c(-1, 1, -1, -1)
#' strategy_multiattribute(ca, cb, v, "TTB")
#' strategy_multiattribute(ca, cb, v, "WADDprob")
#'
#' # multiple item types
#' data(heck2017_raw)
#' strategy_multiattribute(
#'   heck2017_raw[1:10, c("a1", "a2", "a3", "a4")],
#'   heck2017_raw[1:10, c("b1", "b2", "b3", "b4")],
#'   v, "WADDprob"
#' )
#' @export
strategy_multiattribute <- function(cueA, cueB, v, strategy,
                                    c = .50, prior = c(1, 1)) {
  check_cues(cueA, cueB, v)

  if (length(strategy) != 1) {
    # multiple strategies
    strat.list <- lapply(strategy, function(s) {
      strategy_multiattribute(cueA, cueB, v, s, c, prior)
    })
    names(strat.list) <- strategy
    return(strat.list)
  } else if (!is.null(dim(cueA))) {
    # multiple item types
    if (!is.matrix(v)) {
      v <- matrix(v, nrow(cueA), length(v), byrow = TRUE)
    }
    if (length(c) == 1) {
      c <- rep(c, nrow(cueA))
    }
    tmp <- mapply(strategy_multiattribute,
      cueA = as.list(data.frame(t(cueA))),
      cueB = as.list(data.frame(t(cueB))),
      v = as.list(data.frame(t(v))),
      c = c,
      MoreArgs = list(strategy = strategy, prior = prior),
      SIMPLIFY = FALSE
    )
    item.list <- tmp[[1]]
    if (strategy == "baseline") {
      item.list$pattern <- 1:length(tmp)
    } else {
      item.list$pattern <- sapply(tmp, "[[", "pattern")
    }
    names(item.list$pattern) <- rownames(cueA)
    return(item.list)
  } else {
    # single strategy and item type
    pred <- switch(strategy,
      "baseline" = NA,
      "TTBprob" = {
        o <- order(v, decreasing = TRUE)
        diff <- cueB[o] - cueA[o]
        d2 <- diff[diff != 0]
        cnt <- which.max(diff != 0)
        ifelse(length(d2) == 0, 0, 1 / cnt * sign(d2[1]))
      },
      "TTB" = {
        cnt <- strategy_multiattribute(cueA, cueB, v, "TTBprob")$pattern
        sign(cnt)
      },
      "EQW" = {
        cnt_diff <- sum(cueB == 1) - sum(cueA == 1)
        sign(cnt_diff)
      },
      "WADDprob" = {
        logodds <- log(v / (1 - v))
        sum(logodds[cueA == -1 & cueB == 1]) -
          sum(logodds[cueA == 1 & cueB == -1])
      },
      "WADD" = {
        ev <- strategy_multiattribute(cueA, cueB, v, "WADDprob")$pattern
        sign(ev)
      },
      "GUESS" = {
        0.0
      },
      {
        stop("Strategy not supported.")
      }
    )
    pred.list <- list(
      pattern = pred,
      c = ifelse(strategy == "baseline", 1, c),
      ordered = grepl("prob", strategy),
      prior = prior,
      label = strategy
    )
    class(pred.list) <- "strategy"
    return(pred.list)
  }
}

# ' @rdname print
#' @export
print.strategy <- function(x, ...) {
  p <- x$pattern
  lp <- length(x$pattern)
  e_idx <- get_error_idx(x$pattern)
  cat("## Strategy prediction ", x$label, "\n")
  cat("## (length of 'pattern': ", lp,
    "; errors ordered: ", x$ordered, " ; prior = Beta(",
    paste(x$prior, collapse = ","), ")\n",
    sep = ""
  )
  print(head(data.frame(
    Pattern = p,
    Prediction = ifelse(p < 0, "A", ifelse(p > 0, "B", "GUESS")),
    Error = ifelse(e_idx == .5, .5,
      paste0("e", e_idx, " <= ", x$c)
    )
  )))
  if (lp > 6) cat("... [", lp - 6, " predictions omitted]", sep = "")
}


#' Transform Pattern of Predictions to Polytope
#'
#' Transforms ordered item-type predictions to polytope definition.
#' This allows to use Monte-Carlo methods to compute the Bayes factor
#' if the number of item types is large (\code{\link{bf_binom}}).
#'
#' @param strategy a decision strategy returned by \code{\link{strategy_multiattribute}}.
#'
#' @details
#' Note: Only works for models without guessing predictions and
#' without equality constraints (i.e., requires separate error probabilities per item type)
#'
#' @return a list containing the matrix \code{A} and the vector \code{b}
#'      that define a polytope via \code{A*x <= b}.
#' @examples
#' # strategy:  A,B,B,A   e2<e3<e4<e1<.50
#' strat <- list(
#'   pattern = c(-1, 4, 3, -2),
#'   c = .5, ordered = TRUE,
#'   prior = c(1, 1)
#' )
#' pt <- strategy_to_Ab(strat)
#' pt
#'
#' # compare results to encompassing BF method:
#' b <- list(
#'   pattern = 1:4, c = 1,
#'   ordered = FALSE, prior = c(1, 1)
#' )
#' k <- c(2, 20, 18, 0)
#' n <- rep(20, 4)
#' m1 <- strategy_postprob(k, n, list(strat, b))
#' log(m1[1] / m1[2])
#' bf_binom(k, n, pt$A, pt$b, log = TRUE)
#' @export
strategy_to_Ab <- function(strategy) {
  check_strategy(strategy)
  pattern <- strategy$pattern
  c <- strategy$c
  # pu <- rev(get_error_unique(pattern))
  n_error <- get_error_number(pattern)
  if (any(pattern == 0) ||
    length(unique(pattern)) != n_error || n_error < 2) {
    stop("Not working for guesses or if there are less errors than item types.")
  }

  # error in pattern negative: 0<P(b)<c ; positive: c<P(b)<1
  abs_p <- abs(pattern)
  s <- sign(pattern)
  A <- rbind(diag(s), diag(-s))
  colnames(A) <- paste0("p", abs(pattern))
  b <- c(
    ifelse(s == -1, 0, 1),
    ifelse(s == -1, c, -c)
  )

  if (strategy$ordered) {
    # linear order constraints:
    Ao <- matrix(0, factorial(n_error - 1), n_error)
    bo <- rep(0, factorial(n_error - 1))
    cnt <- 0
    for (i in 1:(n_error - 1)) {
      for (j in (i + 1):n_error) {
        cnt <- cnt + 1
        Ai <- pattern[i] < 0
        Aj <- pattern[j] < 0
        ij <- abs_p[i] > abs_p[j] # TRUE <=> ei<ej
        if (Ai && Aj) { # A/A
          Ao[cnt, c(i, j)] <- ifelse(rep(ij, 2), # <=> bi<bj
            c(1, -1),
            c(-1, 1)
          )
        } else if (!Ai && !Aj) { # B/B
          Ao[cnt, c(i, j)] <- ifelse(rep(ij, 2), # <=> 1-bi<1-bj
            c(-1, 1),
            c(1, -1)
          )
        } else if (Ai && !Aj) { # A/B
          Ao[cnt, c(i, j)] <- ifelse(rep(ij, 2),
            c(1, 1), # <=> bi<1-bj
            c(-1, -1)
          ) # <=> bi>1-bj
          bo[cnt] <- ifelse(ij, 1, -1)
        } else if (!Ai && Aj) { # B/A
          Ao[cnt, c(i, j)] <- ifelse(rep(ij, 2),
            c(-1, -1), # <=> 1-bi<bj
            c(1, 1)
          ) # <=> 1-bi>bj
          bo[cnt] <- ifelse(ij, -1, 1)
        }
      }
    }
    A <- rbind(A, Ao)
    b <- c(b, bo)
  }
  list(A = A, b = b)
}
