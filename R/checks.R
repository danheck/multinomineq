
check_bnp <- function(b = NULL, n = NULL, prediction = NULL){
  if(any(n<0) || any(n != round(n)))
    stop("'n' must contain positive integers.")

  if(any(b<0) || any(b>n) || any(b != round(b)))
    stop("'b' must contain positive integers.")

  if(any(is.na(prediction))) # || any(prediction<0) || any(prediction>1) )
    stop("'prediction' must not contain missings.")


  l <- length(n)
  if( l != length(b) || l != length(b) || l != length(prediction))
    stop("Length of 'b', 'n', and 'prediction' does not match.")
}


check_cnml <- function(cnml, n = NULL){
  if (is.list(cnml) && identical(names(cnml), c("cnml", "n", "prediction", "c", "luck", "time")) ){
    if(!missing(n) && !identical(n, cnml$n))
      stop("sample size 'n' not identical for 'cnml'")

    check_luck(cnml$luck)
  } else if (is.list(cnml)) {
    sapply(cnml, check_cnml)
  } else {
    stop("Structure of 'cnml' does not fit to that returned by compute_cnml().")
  }
}

check_luck <- function (luck){
  if (length(luck) != 2 || any(luck <= 0))
    stop ("cnml$luck must be a vector with two positive values")
}
