############ CONSTANTS
# to avoid numerical issues: lower value for probability parameters
BOUND <- 1e-10
MIN_LL <- - 1e300



drop_fixed <- function(k, options = rep(2, length(k) / 2)){
  idx <- cumsum(options)
  if (is.null(dim(k)) || length(dim(k)) == 1){
    k[- idx]
  } else {
    k[,- idx, drop = FALSE]
  }
}

add_fixed <- function(k, sum = 1, options = rep(2, length(k))){
  if (length(sum) == 1)
    sum <- rep(sum, length(k))
  oo <- rep(1:length(options), options - 1)
  k_all <- c()
  for(i in seq_along(options)){
    ko <- k[oo == i]
    k_all <- c(k_all, ko, sum[i] - sum(ko))
  }
  names(k_all) <- index_mult(options)
  k_all
}

k_to_prob <- function(k, options = rep(2, length(k) / 2)){
  oo <- rep(1:length(options), options)
  n <- tapply(k, oo, sum)[oo]
  ml <- k / n
  # ml[is.na(ml)] <- runif(sum(is.na(ml)))
  drop_fixed(ml, options)
}

index_bin <- function(k){
  I <- length(k)
  paste0("p", rep(1:I, each = 2), 1:2)
}

index_mult <- function(options){
  I <- length(options)
  j <- unlist(sapply(options, function(x) seq(1,x)))
  paste0("p", rep(1:I, options), j)
}
