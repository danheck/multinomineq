
check_kn <- function(k = NULL, n = NULL){
  stopifnot(all(n >= 0))
  if (!is.null(dim(k)) && length(dim(k)) == 2)
    stopifnot(ncol(k) == length(n), all(k >= 0), all(t(k) <= n))
  else
    stopifnot(length(k) == length(n), all(k >= 0), all(k <= n))

}


check_start <- function(start, A, b, interior = FALSE){
  if (length(start) != ncol(A))
    stop ("'start' must have the same length as the number of columns of 'A'.")

  if (!interior && !all(A %*% start <= b))
    stop ("'start' must be in the convex polytope:  A*start <= b")

  if (interior && !all(A %*% start < b))
    stop ("'start' must be in the interior (not at the boundaries) of convex polytope:  A*start < b")
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

check_prob <- function  (prob){
  if (is.null(dim(prob)))
    stop("'prob' must be a matrix (usually the posterior samples).")
  if (any(prob < 0, prob > 1))
    stop("'prob' must contain probabilities in [0,1].")
}
check_probkn <- function (prob, k, n){
  check_kn(k, n)
  check_prob(prob)
  if (ncol(prob) != length(k))
    stop("length(k)  must be identical to  ncol(prob).")
}

check_probko <- function (prob, k, options, p_drop = TRUE){
  check_ko(k, options)
  check_prob(prob)
  if (p_drop && ncol(prob) != sum(options - 1))
    stop("If p_drop==TRUE: Number of probabilities must be identical to  sum(options-1).")
  if (!p_drop && ncol(prob) != sum(options))
    stop("If p_drop==FALSE: Number of probabilities must be identical to  sum(options).")
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
  if (!missing(prior) && !is.null(prior) &&
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

check_ko <- function(k, options, label = "k"){
  check_o(options)
  if (length(k) != sum(options))
    stop("'", label, "' must have the same length as sum(options).")
}

check_Abokprior <- function (A, b, options, k, prior = rep(1, length(k))){
  check_Ab(A, b, options)
  check_ko(k, options)
  if (is.null(prior) || length(prior) != length(k) ||
      !is.numeric(prior) || any(prior <0))
    stop ("'prior' must contain nonnegative numbers and have the same length as 'k'.")
}



check_Ab <- function(A, b, options = rep(2, ncol(A))){
  if (length(options) == 1)
    options <- rep(options, ncol(A) / (options - 1))
  if (is.null(dim(A)) ||  nrow(A) != length(b))
    stop("'A' must be a matrix with number of rows equal to the length of 'b'.")
  if (!is.numeric(A) || !is.numeric(b) || any(is.na(A)) || any(is.na(b)))
    stop ("'A' and 'b' must be numeric.")
  check_ko(rep(0, sum(options)), options)
  if (ncol(A) != sum(options - 1))
    stop ("The number of columns in 'A' must be identical to sum(options-1).' ")
  options
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

check_V <- function(V, options = 2){
  if(is.null(dim(V)) ||  any(V < 0, V > 1))
    stop("The vertex representation 'V' must be provided as a numeric matrix with values in [0,1].")

  if (nrow(unique(V)) != nrow(V))
    warning("The matrix 'V' contains identical vertices (rows).\n",
            "  This makes estimation and testing functions unstable and less efficient.\n",
            "  Please remove redundant vertices, e.g., by using:  unique(V)")

  # if (length(options) == 1)
  #   options <- rep(options, ncol(V) / (options - 1))
  # else if (sum(options - 1) != ncol(V))
  #   stop ("V and options do not match: sum(options - 1) != ncol(V) ")
  rep_options(options, V[1,])
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
  if (missing(steps) || is.null(steps))
    steps <- 1:nrow(A)
  if (any(steps <= 0) || any(steps > nrow(A)) || any(steps != round(steps)))
    stop("'steps' must be a vector with positive integers not larger",
         "\n  than the number of rows of the matrix 'A'.")
  unique(sort(c(steps, nrow(A))))
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

check_Mminmax <- function(M, cmin = 0, maxiter = 100, steps){
  stopifnot(all(M > 0), all(M == round(M)))
  if(!missing(cmin) && length(cmin) != 1 || cmin < 0 ||  cmin != round(cmin))
    stop("'cmin' must be a nonnegative integer.")
  if(!missing(maxiter) &&
     length(maxiter) != 1 || maxiter < 1 ||  maxiter != round(maxiter))
    stop("'maxiter' must be a positive integer.")
  if (!missing(steps) && !length(M) %in% c(1, 2, length(steps) + 1:2))
    stop("'M' must be of length 1, 2, or, length(steps). \n",
         "  (identical number of iterations for all steps)")
}

check_count <- function(count){
  if (is.null(dim(count)))
    count <- matrix(count, 1, dimnames = list(NULL, names(count)))

  # if(!is.numeric(count) || ! all(c("count", "M") %in% colnames(count)))
  #   stop ("'prior' / 'posterior' must be a matrix with the column names 'count' and 'M")
  cc <- count[,c("count","M")]
  if(any(cc < 0) ||  any(cc != round(cc)))
    stop ("'count' and must contain positive integers.")
  count
}

check_io <- function(inside, options){
  x <- c(rpdirichlet(1, rep(1.5, sum(options)), options))
  if(is.function(inside)){
    tryCatch (inside_output <- inside(x),
              error = function(e)
                stop("The function 'inside' must be valid for vector of length ", sum(options - 1)))
    stopifnot(inside_output %in% c(0, 1))
  } else if(class(inside) == "XPtr"){
    RcppXPtrUtils::checkXPtr(inside, type = "SEXP", args = "NumericVector")
    tryCatch (inside_output <- call_xptr(inside, x),
              error = function(e)
                stop("The C++ function 'inside' defined via RcppXPtrUtils::call_xptr",
                     "\n  must be valid for vector of length ", sum(options - 1)))
  } else {
    stop("'inside' must be an R function or a C++ pointer to a function\n",
         "generated via RcppXPtrUtils::cppXPtr(code) .")
  }
  if (!inside_output %in% c(0,1))
    stop ("The function 'inside' should return values 0/1  or TRUE/FALSE.")
}

