
index_bin <- function(k){
  I <- length(k)
  paste0("i", rep(1:I, each = 2), 1:2)
}

index_mult <- function(options){
  I <- length(options)
  j <- unlist(sapply(options, function(x) seq(1,x)))
  paste0("i", rep(1:I, options), j)
}
