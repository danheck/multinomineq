#' Posterior Predictive p-Values
#'
#' Uses posterior samples to get posterior-predicted frequencies and compare
#' the Pearson's X^2 statistic for (1) the observed frequencies vs. (2) the posterior-predicted frequencies.
#'
#' @inheritParams count_binomial
#' @inheritParams count_multinomial
#' @param theta posterior samples of the choice probabilities.
#'     See \code{\link{sampling_binomial}} and \code{\link{sampling_multinomial}}
#' @param M number of subsamples from \code{theta}
#' @seealso \code{\link{sampling_binomial}} and \code{\link{sampling_multinomial}} to get posterior samples
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
#' ppp_multinomial(theta, c(1,0,9), 3)
#' @export
ppp_binomial <- function(theta, k, n, M = 2000){
  check_thetakn(theta, k, n)
  M <- min(nrow(theta), M)
  tt <- theta[sample(nrow(theta), M),]
  n_mat <- matrix(n, M, length(n), byrow = TRUE)
  k_obs <- matrix(k, M, length(k), byrow = TRUE)
  k_pred <- matrix(rbinom(M*length(n), n, tt), M, byrow = TRUE)

  X2_pred <- X2(e = cbind(tt, 1 - tt) * n,
                o = cbind(k_pred, n_mat - k_pred))
  X2_obs <- X2(e = cbind(tt, 1 - tt) * n,
               o = cbind(k_obs, n_mat - k_obs))
  c("X2_obs" = mean(X2_obs), "X2_pred" = mean(X2_pred),
    "ppp" = mean(X2_obs < X2_pred))
}

#' @rdname ppp_binomial
#' @export
ppp_multinomial <- function(theta, k, options, M = 2000){
  check_thetako(theta, k, options)
  M <- min(nrow(theta), M)
  tf <- free_to_full(theta[sample(nrow(theta), M),], options)
  n <- tapply(k, rep(seq_along(options), options), sum)
  n_mat <- matrix(rep(n, options), M, sum(options), byrow= TRUE)
  k_obs <- matrix(k, M, length(k), byrow = TRUE)
  k_pred <- pred_mult(tf, k, options)

  X2_pred <- X2(e = tf * n_mat, o = k_pred)
  X2_obs <- X2(e = tf * n_mat, o = k_obs)
  c("X2_obs" = mean(X2_obs), "X2_pred" = mean(X2_pred),
    "ppp" = mean(X2_obs < X2_pred))
}

# Pearson's X^2
X2 <- function(e, o){
  if (is.null(dim(e))){
    sum( (o - e)^2 / e )
  } else {
    rowSums( (o - e)^2 / e )
  }
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

