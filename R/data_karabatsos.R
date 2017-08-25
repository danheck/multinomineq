#' Posterior Predictive p-Values
#'
#' Uses posterior samples to get posterior-predicted frequencies and compare
#' the Pearson's X^2 statistic for (1) the observed frequencies vs. (2) the posterior-predicted frequencies.
#'
#' @inheritParams rpmultinom
#' @inheritParams count_binomial
#' @inheritParams count_multinomial
#' @param by optional: a vector of the same length as \code{k} indicating factor levels
#'     by which the posterior-predictive checks should be split (e.g., by item sets).
#' @param M number of subsamples from \code{theta}
#' @seealso \code{\link{sampling_binomial}}/\code{\link{sampling_multinomial}} to get posterior samples and \code{\link{rpbinom}}/\code{\link{rpmultinom}} to get posterior-predictive samples.
#' @template ref_myung2005
#'
#' @examples
#' # uniform samples:  p<.10
#' theta <- matrix(runif(300*3, 0, .1), 300)
#' n <- rep(10, 3)
#' ppp_binomial(theta, c(1,2,0), n)  # ok
#' ppp_binomial(theta, c(5,4,3), n)  # misfit
#'
#' # multinomial (ternary choice)
#' theta <- matrix(runif(300*2, 0, .05), 300)
#' ppp_multinomial(theta, c(1,0,9), 3)  # ok
#' @export
ppp_binomial <- function(theta, k, n, by){
  if (nrow(theta) < M) {
    warning("Not enough samples, 'M' is adjusted.")
    M <- nrow(theta)
  }
  if (length(n) == 1)
    n <- rep(n, length(k))
  check_thetakn(theta, k, n)
  if (!missing(by)){
    if (length(by) != length(k))
      stop("length of 'by' does not match length of 'k'.")
    levels <- unique(by)
    res <- matrix(NA, length(levels), 3,
                  dimnames = list(levels, c("X2_obs", "X2_pred", "ppp")))
    for (i in 1:length(levels)){
      sel <- by == levels[i]
      res[i,] <- ppp_binomial(theta[,sel, drop = FALSE], k[sel], n[sel], M = M)
    }

  } else {
    tt <- theta[sample(nrow(theta), M),,drop = FALSE]
    res <- ppp_bin(theta[sample(nrow(theta), M),,drop = FALSE], k, n)
  }
  res
}

#' @rdname ppp_binomial
#' @export
ppp_multinomial <- function(theta, k, options, by, M = 2000){
  check_thetako(theta, k, options)
  if (nrow(theta) < M) {
    warning("Not enough samples, 'M' is adjusted.")
    M <- nrow(theta)
  }
  tf <- free_to_full(theta[sample(nrow(theta), M),], options)
  n <- tapply(k, rep(seq_along(options), options), sum)
  n_mat <- matrix(rep(n, options), M, sum(options), byrow= TRUE)
  k_obs <- matrix(k, M, length(k), byrow = TRUE)
  k_pred <- rpmultinom(tf, get_n(k, options), options)

  X2_pred <- X2(e = tf * n_mat, o = k_pred)
  X2_obs  <- X2(e = tf * n_mat, o = k_obs)
  c("X2_obs" = mean(X2_obs),
    "X2_pred" = mean(X2_pred),
    "ppp" = mean(X2_obs <= X2_pred))
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

free_to_full <- function(theta, options){
  oo <- rep(1:length(options), options - 1)
  if (is.null(dim(theta))){
    p_J <- 1 - tapply(theta, oo, function(x) sum(x))
    full <- unlist(mapply(c, by(theta, oo, list), p_J))
  } else {
    p_J <- by(data.frame(t(theta)), oo, function(x) 1 - colSums(x))
    o_list <- mapply(rbind, by(t(theta), oo, list), p_J, SIMPLIFY = FALSE)
    full <- t(do.call("rbind", o_list))
    colnames(full) <- index_mult(options)
  }
  full
}

# bin_to_mult <- function(k){
#   list("k" = rep(k, each = 2), "options" = rep(2, length(k)))
# }

