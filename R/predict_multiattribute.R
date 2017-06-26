
#' Get Predictions of Common Strategies
#'
#' Returns a vector of predictions for well-known choice strategies (TTB, WADD)
#' @param cueA cue values of Option A (-1/+1 = negative/positive; 0 = missing).
#'        If a matrix is provided, each row defines one item type.
#' @param cueB cue values of Option B (see \code{cueA}).
#' @param v cue validities: probabilities that cues lead to correct decision.
#'        Must be of the same length as the number of cues.
#' @param strategy strategy label, e.g., \code{"TTB"}, \code{"WADD"}, or \code{"WADDprob"}. See details.
#' @return Predidcted choice pattern are encoded as negative = Option A, positive = Option B, and 0 = guessing.
#'         Different integers are used if error probabilities are assumed to be different.
#' @examples
#' ca <- c(1, -1, -1, 1)
#' cb <- c(-1, 1, -1, -1)
#' v <- c(.9, .8, .7, .6)
#' predict_multiattribute(ca, cb, v, "TTBprob")
#' predict_multiattribute(ca, cb, v, "TTB")
#' predict_multiattribute(ca, cb, v, "WADDprob")
#' predict_multiattribute(ca, cb, v, "WADD")
#' predict_multiattribute(ca, cb, v, "EQW")
#' predict_multiattribute(ca, cb, v, "GUESS")
#' @export
predict_multiattribute <- function (cueA, cueB, v, strategy){
  UseMethod("predict_multiattribute", cueA)
}

#' @rdname predict_multiattribute
#' @export
predict_multiattribute.default <- function (cueA, cueB, v, strategy){
  check_cues(cueA, cueB, v)
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
  if (strategy %in% c("WADDprob", "TTBprob"))
    attr(pred, "ordered") <- TRUE
  else
    attr(pred, "ordered") <- NULL
  pred
}

#' @rdname predict_multiattribute
#' @export
predict_multiattribute.matrix <- function (cueA, cueB, v, strategy){
  check_cues(cueA, cueB, v)
  if (!is.matrix(v)) v <- matrix(v, nrow(cueA), length(v), byrow = TRUE)
  res <- mapply(predict_multiattribute,
                cueA = as.list(data.frame(t(cueA))),
                cueB = as.list(data.frame(t(cueB))),
                v = as.list(data.frame(t(v))),
                MoreArgs = list(strategy = strategy))
  names(res) <- rownames(cueA)
  res
}

#' @rdname predict_multiattribute
#' @export
predict_multiattribute.data.frame <- function (cueA, cueB, v, strategy){
  cueA <- as.matrix(cueA)
  predict_multiattribute(cueA, cueB, v, strategy)
}
