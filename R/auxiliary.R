############ CONSTANTS
# to avoid numerical issues: lower value for probability parameters
BOUND <- 1e-10
MIN_LL <- - 1e100

# reduce memory load by counting in batches:
BATCH <- 5000


############################# S3 methods: drop / add fixed dimensions
#' Converts Binary to Multinomial Frequencies
#'
#' Converts the number of "hits" in the binary choice format to the observed
#' frequencies across  for all response categories (i.e., the multinomial format).
#'
#' @inheritParams count_binom
#'
#' @details
#' In \code{multinomineq}, binary choice frequencies are represented by the number
#' of "hits" for each item type/condition (the vector \code{k}) and by the total
#' number of responses per item type/condition (the scalar or vector \code{n}).
#'
#' In the multinomial format, the vector \code{k} includes all response categories
#' (not only the number of "hits"). This requires to define a vector \code{options},
#' which indicates how many categories belong to one item type/condition (since
#' the total number of responses per item type is fixed).
#'
#' @examples
#' k <- c(1, 5, 8, 10)
#' n <- 10
#' binom_to_multinom(k, n)
#'
#' @export
binom_to_multinom <- function(k, n){
  if (length(n) == 1) n <- rep(n, length(k))
  check_kn(k, n)
  options <- rep(2, length(k))
  k_fixed <- add_fixed(k, options, n)
  list(k = k_fixed, options = options)
}

############################# S3 methods: drop / add fixed dimensions

#' Drop or Add Fixed Dimensions for Multinomial Probabilities/Frequencies
#'
#' Switches between two representation of polytopes for multinomial probabilities
#' (whether the fixed parameters are included).
#'
#' @param x a vector (typically \code{k}, \code{n}, or \code{prior}) or
#'   a matrix (typically \code{A} or \code{V}), in which case the fixed dimensions
#'   are dropped/added column-wise.
#' @inheritParams count_multinom
#'
#' @examples
#' ######## bi- and trinomial (a1,a2, b1,b2,b3)
#' # vectors with frequencies:
#' drop_fixed(c(3,7,  4,1,5), options = c(2,3))
#' add_fixed (c(3,    4,1),   options = c(2,3),
#'            sum = c(10, 10))
#'
#' # matrices with probabilities:
#' V <- matrix(c(1,  0,  0,
#'               1, .5, .5,
#'               0,  1,  0), 3, byrow = TRUE)
#' V2 <- add_fixed(V, options = c(2,3))
#' V2
#' drop_fixed(V2, c(2,3))
#' @export
drop_fixed <- function(x, options = 2)
  UseMethod("drop_fixed", x)

#' @export
drop_fixed.matrix <- function(x, options = 2){
  options <- rep_options(options, x[1,], p_drop = FALSE)
  xn <- x[,- cumsum(options), drop = FALSE]
  if (is.null(colnames(xn)))
    colnames(xn) <- index_mult(options, fixed = FALSE)
  xn
}

#' @export
drop_fixed.data.frame <- function(x, options = 2){
  drop_fixed.matrix(x, options)
}

#' @export
drop_fixed.default <- function(x, options = 2){
  options <- rep_options(options, x, p_drop = FALSE)
  check_ko(k = x, options = options, label = "x")
  xn <- x[- cumsum(options)]
  if (is.null(names(xn)))
    names(xn) <- index_mult(options, fixed = FALSE)
  xn
}


#' @param sum a vector that determines the fixed sum in each multinomial condition.
#'   By default, probabilities are assumed that sum to one.
#'   If frequencies \code{n} are provided, use \code{sum=n}.
#' @rdname drop_fixed
#' @export
add_fixed <- function(x, options = 2, sum = 1){
  UseMethod("add_fixed", x)
}

# ' @rdname drop_fixed
#' @importFrom stats start end
#' @import coda
#' @export
add_fixed.mcmc <- function(x, options = 2, sum = 1){
  as.mcmc(add_fixed.matrix(x, options = options, sum = sum),
          start = start(x), thin = thin(x), end = end(x))
}

# ' @rdname drop_fixed
#' @export
add_fixed.matrix <- function(x, options = 2, sum = 1){
  options <- rep_options(options, x[1,])
  if (length(sum) == 1)
    sum <- rep(sum, length(options))
  xn <- matrix(NA, nrow(x), sum(options))
  for (i in seq_along(options)){
    k <- sum(options[0:(i-1)])       + seq(1, options[i] - 1)
    l <- sum((options - 1)[0:(i-1)]) + seq(1, options[i] - 1)
    xn[,k] <- x[,l]
    rs <- rowSums(xn[,k, drop = FALSE])
    message <- paste0("Within multinomial condition ", i, ", row sums of matrix must be: ",
                      " <= ",sum[i],"  and  >= 0 ")
    if (any(rs < 0 | rs > sum[i]))
      stop(message)
    xn[,max(k) + 1] <- sum[i] - rs
  }
  colnames(xn) <- index_mult(options, fixed = TRUE)
  xn
}

# ' @rdname drop_fixed
#' @export
add_fixed.data.frame <- function(x, options = 2, sum = 1){
  add_fixed.matrix(x, options, sum)
}

# ' @rdname drop_fixed
#' @export
add_fixed.default <- function(x, options = 2, sum = 1){
  add_fixed.matrix(t(x), options, sum)[1,]
  # check_ko(k, options)
  # oo <- rep(1:length(options), options - 1)
  # xn <- c()
  # for(i in seq_along(options)){
  #   ko <- x[oo == i]
  #   k_all <- c(k_all, ko, sum[i] - sum(ko))
  # }
  # names(k_all) <- index_mult(options)
  # k_all
}



index_free_to_fixed <- function(i, options){
  option <- rep(1:length(options), options - 1)[i]
  i + option - 1
}

rep_options <- function(options, x, p_drop = TRUE){
  if (length(options) == 1){
    times <- length(x) / (options - p_drop)
    if (times != round(times))
      stop("Check input: The length/number of columns of 'x' is not a multiple of (options-1)")
    options <- rep(options, times)
  } else if (length(x) != sum(options - p_drop)){
    stop("Length of 'options' does not match number of frequencies/probabilities/columns.")
  }
  options
}

k_to_prob <- function(k, options = rep(2, length(k) / 2)){
  oo <- rep(1:length(options), options)
  n <- tapply(k, oo, sum)[oo]
  ml <- k / n
  # ml[is.na(ml)] <- runif(sum(is.na(ml)))
  drop_fixed(ml, options)
}

# meaningful labels for column indices of A and V / names of k/n
index_bin <- function(k){
  I <- length(k)
  paste0("p", rep(1:I, each = 2), "_", 1:2)
}
index_mult <- function(options, fixed = TRUE){
  I <- length(options)
  if (fixed){
    j <- unlist(sapply(options, function(x) seq(1,x)))
    labels <- paste0("p", rep(1:I, options), "_", j)
  } else {
    j <- unlist(sapply(options - 1, function(x) seq(1,x)))
    labels <- paste0("p", rep(1:I, options - 1), "_", j)
  }
  labels
}
