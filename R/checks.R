
check_bnp <- function(k = NULL, n = NULL, prediction = NULL){
  if(any(n < 0) || any(n != round(n)))
    stop("'n' must contain positive integers.")

  if(any(k < 0) || any(k > n) || any(k != round(k)))
    stop("'k' must contain positive integers smaller or equal to 'n'.")

  if(any(is.na(prediction))) # || any(prediction<0) || any(prediction>1) )
    stop("'prediction' must not contain missings.")


  l <- length(n)
  if( l != length(k) || l != length(k) || l != length(prediction))
    stop("Length of 'k', 'n', and 'prediction' does not match.")
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
    stop ("cnml$luck must be a vector with two positive values.")
}

check_cues <- function (cueA, cueB, v){
  if (is.matrix(cueA) || is.data.frame(cueA)){
    if (!identical(dim(cueA), dim(cueB)))
      stop("Size of matrices for 'cueA' and 'cueB' must match.")
  } else {
    if (length(cueA) != length(cueB) || length(cueA) != length(v))
      stop("Length of 'cueA', 'cueB', and 'v' must match.")
  }
  if (!all(unlist(c(cueA,cueB)) %in% c(-1,0,1)))
    stop ("Cues must have values -1/0/+1.")
  if (any(v < 0, v > 1))
    stop ("Validities 'v' must be in [0,1].")
}

check_knpcp <- function(k, n, pred, c, prior = c(1, 1)){
  check_bnp(k, n, pred)
  if (!is.numeric(c) || length(c) != 1 || c < 0 || c > 1)
    stop("'c' must be in the interval [0,1].")
  if (!missing(prior) &&
      (length(prior) != 2 || any(prior < 0) || !is.numeric(prior)))
    stop("'prior' must be a numeric vector with two postive numbers.")
}
