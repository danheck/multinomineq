#' Posterior Predictive p-Values
#'
#' Uses posterior samples to get posterior-predicted frequencies and compare
#' the Pearson's X^2 statistic for (1) the observed frequencies vs. (2) the posterior-predicted frequencies.
#'
#' @inheritParams rpmultinom
#' @inheritParams count_binom
#' @inheritParams count_multinom
#' @param by optional: a vector of the same length as \code{k} indicating factor levels
#'     by which the posterior-predictive checks should be split (e.g., by item sets).
# ' @param M number of subsamples from \code{prob}
#' @seealso \code{\link{sampling_binom}}/\code{\link{sampling_multinom}} to get posterior samples and \code{\link{rpbinom}}/\code{\link{rpmultinom}} to get posterior-predictive samples.
#' @template ref_myung2005
#'
#' @examples
#' # uniform samples:  p<.10
#' prob <- matrix(runif(300*3, 0, .1), 300)
#' n <- rep(10, 3)
#' ppp_binom(prob, c(1,2,0), n)  # ok
#' ppp_binom(prob, c(5,4,3), n)  # misfit
#'
#' # multinomial (ternary choice)
#' prob <- matrix(runif(300*2, 0, .05), 300)
#' ppp_multinom(prob, c(1,0,9), 3)  # ok
#' @export
ppp_binom <- function(prob, k, n, by){
  if (length(n) == 1)
    n <- rep(n, length(k))
  check_probkn(prob, k, n)
  if (!missing(by)){
    if (length(by) != length(k))
      stop("length of 'by' does not match length of 'k'.")
    levels <- unique(by)
    res <- matrix(NA, length(levels), 3,
                  dimnames = list(levels, c("X2_obs", "X2_pred", "ppp")))
    for (i in 1:length(levels)){
      sel <- by == levels[i]
      res[i,] <- ppp_binom(prob[,sel, drop = FALSE], k[sel], n[sel])
    }

  } else {
    res <- ppp_bin(matrix(prob, ncol = length(k)), k, n)
  }
  res
}

#' @rdname ppp_binom
#' @export
ppp_multinom <- function(prob, k, options, p_drop = TRUE){
  check_probko(prob, k, options, p_drop = p_drop)
  if (p_drop)
    prob <- add_fixed(prob, options)
  ppp_mult(prob, k, options)

  ### TODO: by factor (necessary within C++!  => multinomial sampling)
  # idea: provide integer vector and use   find(by == i)  in RcppArmadillo
}

# Pearson's X^2
X2 <- function(e, o){
  if (is.null(dim(e))){
    sum( (o - e)^2 / e )
  } else {
    rowSums( (o - e)^2 / e )
  }
}

get_n <- function(k, options){
  oo <- rep(1:length(options), options)
  tapply(k, oo, sum)
}

pred_mult <- function(tf, k, options){
  oo <- rep(1:length(options), options)
  n <- tapply(k, oo, sum)
  k_pred <- matrix(NA, nrow(tf), sum(options))
  for(i in 1:nrow(tf)){
    k_pred[i,] <- unlist(mapply(rmultinom, size = n, prob = split(tf[i,], oo),
                                MoreArgs = list(n = 1)))
  }
  k_pred
}

# redundant due to add_fixed!
# free_to_full <- function(prob, options){
#   oo <- rep(1:length(options), options - 1)
#   if (is.null(dim(prob))){
#     p_J <- 1 - tapply(prob, oo, function(x) sum(x))
#     full <- unlist(mapply(c, by(prob, oo, list), p_J))
#   } else {
#     p_J <- by(data.frame(t(prob)), oo, function(x) 1 - colSums(x))
#     o_list <- mapply(rbind, by(t(prob), oo, list), p_J, SIMPLIFY = FALSE)
#     full <- t(do.call("rbind", o_list))
#     colnames(full) <- index_mult(options)
#   }
#   full
# }

# bin_to_mult <- function(k){
#   list("k" = rep(k, each = 2), "options" = rep(2, length(k)))
# }

