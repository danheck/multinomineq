
check_kn <- function(k = NULL, n = NULL){
  if(length(k) != length(n))
    stop("Length of 'k' and 'n' does not match.")

  if(any(n < 0) || any(n != round(n)))
    stop("'n' must contain positive integers.")

  if(any(k < 0) || any(k > n) || any(k != round(k)))
    stop("'k' must contain positive integers smaller or equal to 'n'.")
}


check_prior <- function (prior){
  if (is.null(prior) || length(prior) != 2 ||
      !is.numeric(prior) || any(prior <= 0))
    stop ("prior must be a vector with two positive values.")
}

check_knp <- function(k = NULL, n = NULL, pattern = NULL){
  check_kn(k, n)
  if(any(is.na(pattern)))
    stop("'pattern' must not contain missings.")

  l <- length(n)
  if( l != length(pattern))
    stop("Length of 'k'/'n', and 'pattern' does not match.")
}

check_theta <- function  (theta){
  if (is.null(dim(theta)))
    stop("'theta' must be a matrix with posterior samples")
  if (any(theta < 0, theta > 1))
    stop("'theta' must contain probabilities in [0,1].")
}
check_thetakn <- function (theta, k, n){
  check_kn(k, n)
  check_theta(theta)
  if (ncol(theta) != length(k))
    stop("length(k)  must be identical to  ncol(theta).")
}

check_thetako <- function (theta, k, options){
  check_ko(k, options)
  check_theta(theta)
  if (ncol(theta) != sum(options - 1))
    stop("ncol(theta)  must be identical to  sum(options-1).")
}


check_cnml <- function(strategy, n){
  if(is.null(strategy$n) || !identical(n, strategy$n))
    stop("sample size 'n' missing in 'strategy' object or not identical")
  if (is.null(strategy$cnml))
    stop("the NML complexity term needs to be pre-computed with ?compute_cnml")

  check_prior(strategy$prior)
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

check_Abknprior <- function (A, b, k, n, prior = c(1, 1)){
  check_prior(prior)
  check_kn(k, n)
  check_Ab(A, b)
  if (ncol(A) != length(k))
    stop("'A' must be a matrix with number of columns equal to the length of 'k'.")
}

check_o <- function(options)
  if (any(options != round(options), options <0))
    stop ("'options' must contain positive integers.")

check_ko <- function(k, options){
  check_o(options)
  if (length(k) != sum(options))
    stop("'k' must have the same length as sum(options).")
}

check_Abokprior <- function (A, b, options, k, prior = rep(1, length(k))){
  check_Ab(A, b, options)
  check_ko(k, options)
  if (length(prior) != length(k))
    stop("'prior' must have the same length as 'k'.")
  if (any(prior != round(prior), prior <0))
    stop ("'prior' must contain positive integers.")
}



check_Ab <- function(A, b, options = rep(2, ncol(A))){
  if (is.null(dim(A)) ||  nrow(A) != length(b))
    stop("'A' must be a matrix with number of rows equal to the length of 'b'.")
  if (!is.numeric(A) || !is.numeric(b) || any(is.na(A)) || any(is.na(b)))
    stop ("'A' and 'b' must be numeric.")
  check_ko(rep(0, sum(options)), options)
  if (ncol(A) != sum(options - 1))
    stop ("The number of columns in 'A' must be identical to sum(options-1).' ")
}

check_Abx <- function (A, b, x){
  check_Ab(A, b)
  if (is.null(dim(x)) || length(dim(x)) == 1){
    if (length(x) != ncol(A))
      stop ("Probability vector 'x' must have the same length as ncol(A).")
  } else {
    if (ncol(x) != ncol(A))
      stop ("Matrix with vertices 'x' must have the same number of columns as ncol(A).")
  }
}

check_V <- function(V){
  if(is.null(dim(V)) ||  any(V < 0, V > 1))
    stop("The vertex representation 'V' must be provided as a numeric matrix with values in [0,1].")
}

check_Vx <- function (V, x){
  check_V(V)
  if (is.null(dim(x)) || length(dim(x)) == 1){
    if (length(x) != ncol(V))
      stop ("Probability vector 'x' must have the same length as ncol(V).")
  } else {
    if (ncol(x) != ncol(V))
      stop ("Matrix with vertices 'x' must have the same number of columns as ncol(V).")
  }
}

check_stepsA <- function(steps, A){
  if (any(steps <= 0) || any(steps >= nrow(A)) || any(steps != round(steps)))
    stop("'steps' must be a vector with positive integers smaller",
         "\n  than the number of rows of the matrix 'A'.")
}

check_strategy <- function (strategy){
  if (!is.list(strategy))
    stop("'strategy' must be a list (?predict_multiattribute) such as: \n",
         "    list(pattern=c(-1,2,0), c=.5, ordered=TRUE, prior=c(1,1))")
  if (!all(c("pattern", "c","ordered","prior") %in% names(strategy)))
    stop("the list 'strategy' must include named elements, e.g.: \n",
         "    list(pattern=c(-1,2,0), c=.5, ordered=TRUE, prior=c(1,1))")
  check_prior(strategy$prior)
  c <- strategy$c
  if (length(c) != 1 || c < 0 || c > 1)
    stop("'c' must be a single numeric value in [0,1]")
  ordered <- strategy$ordered
  if (is.null(ordered) || length(ordered) != 1 || !is.logical(ordered))
    stop("'strategy$ordered' must be a logical value")
}

check_data_strategy <- function (k, n, strategy){
  check_strategy(strategy)
  check_knp(k, n, strategy$pattern)
}

check_Mbatch <- function(M, batch){
  if(any(M < 0) ||  any(M != round(M)))
    stop("'M' must contain positive integers.")
  if(length(batch) != 1 || batch < 0 ||  batch != round(batch))
    stop("'batch' must be a positive integer.")
}

check_count <- function(count){
  if (!is.list(count) || ! c("count", "M") %in% names(count))
    stop ("'prior' / 'posterior' msut be a list with the entries 'count' and 'M")
  if(any(count$M < 0) ||  any(count$M != round(count$M)))
    stop ("'M' must contain positive integers.")
  if(any(count$count < 0) ||  any(count$count != round(count$count)))
    stop ("'count' must contain positive integers.")
  if (length(count$count) != length(count$M))
    stop ("in 'prior'/'posterior': 'count' and 'M' must have equal length.")
}

check_io <- function(inside, options){
  if (!is.function(inside))
    stop("'inside' must be a function.")
  tryCatch (inside(runif(sum(options))),
            error = function(e)
              stop("The function 'inside' should work for vector of length ", sum(options)))
  if (!inside(runif(sum(options))) %in% c(0,1))
    stop ("The function 'inside' should return 0/1  or TRUE/FALSE.")
}

