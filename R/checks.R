
check_knp <- function(k = NULL, n = NULL, prediction = NULL){
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


check_cnml <- function(strategy, n){
  if(is.null(strategy$n) || !identical(n, strategy$n))
    stop("sample size 'n' missing in 'strategy' object or not identical")
  if (is.null(strategy$cnml))
    stop("the NML complexity term needs to be pre-computed with ?compute_cnml")

  check_luck(strategy$prior)
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
  check_knp(k, n, pred)
  if (!is.numeric(c) || length(c) != 1 || c < 0 || c > 1)
    stop("'c' must be in the interval [0,1].")
  if (!missing(prior) &&
      (length(prior) != 2 || any(prior < 0) || !is.numeric(prior)))
    stop("'prior' must be a numeric vector with two postive numbers.")
}

check_kAb <- function (k, A, b){
  if (is.null(dim(A)) || ncol(A) != length(k) || nrow(A) != length(b))
    stop("'A' must be a matrix with number of rows equal to the length of 'b'",
         "\n  and number of columns equal to the length/column number of 'k'.")
  if (!is.numeric(A) || !is.numeric(b) || any(is.na(A)) || any(is.na(b)))
    stop ("'A' and 'b' must be numeric.")
}

check_stepsA <- function(steps, A){
  if (any(steps <= 0) || any(steps >= nrow(A)) || any(steps != round(steps)))
    stop("'steps' must be a vector with positive integers smaller",
         "\n  than the number of rows of the matrix 'A'.")
}

check_strategy <- function (strategy){
  if (!is.list(strategy))
    stop("'strategy' must be a list")
  if (!all(c("pattern", "c","ordered","prior") %in% names(strategy)))
    stop("'strategy' must have the named elements: \n",
         "    pattern, c, ordered, prior, label")

}

check_data_strategy <- function (k, n, strategy){
  check_strategy(strategy)
  check_knp(k, n, strategy$pattern)
}
