#########################
#' Log-Marginal Likelihood for Decision Strategy
#'
#' Computes the logarithm of the integral over the likelihood function weighted
#' by the prior distribution of the error probabilities.
#'
#' @param k observed frequencies of Option B.
#'          Either a vector or a matrix/data frame (one person per row).
#' @param n vector with the number of choices per item type.
#' @param strategy a list that defines the predictions of a strategy, see\code{\link{strategy_multiattribute}}.
#'
#'
#' @examples
#' k <- c(1,11,18)
#' n <- c(20, 20, 20)
#' # pattern: A, A, B with constant error e<.50
#' strat <- list(pattern = c(-1, -1, 1),
#'               c = .5, ordered = FALSE,
#'               prior = c(1,1))
#' m1 <- strategy_marginal(k, n, strat)
#' m1
#'
#' # pattern: A, B, B with ordered error e1<e3<e2<.50
#' strat2 <- list(pattern = c(-1, 3, 2),
#'                c = .5, ordered = TRUE,
#'                prior = c(1,1))
#' m2 <- strategy_marginal(k, n, strat2)
#' m2
#'
#' # Bayes factor: Model 2 vs. Model 1
#' exp(m2 - m1)
#' @export
strategy_marginal <- function (k, n, strategy){
  check_data_strategy(k, n, strategy)
  prior <- strategy$prior
  c <- strategy$c

  n_error <- get_error_number(strategy$pattern)
  error     <- count_errors(k, n,   strategy$pattern,
                            ordered = strategy$ordered, prob = FALSE)
  adherence <- count_errors(k, n, - strategy$pattern,
                            ordered = strategy$ordered, prob = FALSE)
  pm <- 0

  if (!strategy$ordered || n_error == 1){
    ### includes:
    # baseline model: separate e_i in [0,1]
    # modal-choice model: separate errors    e_i in [0,c]
    # deterministic (constant error): one parameter   e in [0,c]
    s1 <- error + prior[1]
    s2 <- adherence + prior[2]
    # integral & posterior normalization:
    pm <- sum(pbeta(c, s1, s2, log.p = TRUE) +
                lbeta(s1, s2) +
                - lbeta(prior[1], prior[2]))    # prior
    # prior normalization due to truncation:
    if (c != 1)
      pm <- pm - n_error * pbeta(c, prior[1],prior[2], log.p = TRUE)

    ### DEPRECATED baseline / deterministic
    # pm <- sum(lbeta(k + prior[1], n - k + prior[2]) - lbeta(prior[1], prior[2]))
    # pm <- pbeta(c, s1, s2, log.p = TRUE) + lbeta(s1, s2) + # prior*lik (unnormalized)
    #   - pbeta(c, prior[1], prior[2], log.p = TRUE) - lbeta(prior[1], prior[2])
    # (normalization of order-constrained prior)

  } else if (n_error <= 6 && n_error != 0){
    e_idx <- get_error_idx(strategy$pattern)
    # probabilistic (linear order constraint; nested unidimensional integration)
    ff <- function(e, idx = 1){
      if (idx == 1){
        e^(error[1] + prior[1] - 1) * (1-e)^(adherence[1] + prior[2] - 1)
      } else if (idx == 2) {
        e^(error[2] + prior[1] - 1) * (1-e)^(adherence[2] + prior[2] - 1) *
          pbeta(e, error[1] + prior[1], adherence[1] + prior[2])*
          beta(error[1] + prior[1], adherence[1] + prior[2])
      } else {
        e^(error[idx] + prior[1] - 1) * (1-e)^(adherence[idx] + prior[2] - 1)*
          sapply(e, function(ee)
            integrate(ff, lower = 0, upper = ee, idx = idx - 1)$value)
      }
    }
    pm <- log(integrate(ff, 0, c, idx = n_error)$value) +
      lfactorial(n_error) - n_error* (pbeta(c, prior[1], prior[2], log.p = TRUE) +
                                    lbeta(prior[1], prior[2]))
  } else if (n_error > 6){
    stop("For more than 6 item types, use ?as_polytope and ?bf_binom.")
  }

  if (any(strategy$pattern == 0)){
    # GUESS
    pm <- pm + sum(n[strategy$pattern == 0]) * log(.50)
  }
  sum(lchoose(n, k)) + pm
}


#' Strategy Classification: Posterior Model Probabilities
#'
#' Posterior model probabilities for multiple strategies (with equal prior model probabilities).
#'
# ' @inheritParams strategy_cnml
# ' @inheritParams strategy_nml
#' @inheritParams strategy_marginal
#' @param strategies list of strategies. See \link{strategy_multiattribute}
#' @param cores number of processing units for parallel computation.
# @param c upper boundary for parameters/error probabilities.
#   A vector of the same length as \code{strategies}
#   can be provided to use a different upper bound per model.
# @param ordered whether error probabilities are assumed to be ordered.
#   Similar as for \code{c}, a vector of the same length as \code{strategies}
#   can be provided to define which models have ordered error probabilities.
#'
#' @seealso \code{\link{strategy_marginal}} and \code{\link{model_weights}}
#' @examples
#' # pattern 1: A, A, B with constant error e<.50
#' strat1 <- list(pattern = c(-1, -1, 1),
#'                c = .5, ordered = FALSE,
#'                prior = c(1,1))
#' # pattern 2: A, B, B with ordered error e1<e3<e2<.50
#' strat2 <- list(pattern = c(-1, 3, 2),
#'                c = .5, ordered = TRUE,
#'                prior = c(1,1))
#' baseline <- list(pattern = 1:3, c = 1, ordered = FALSE,
#'                  prior = c(1,1))
#'
#' # data
#' k <- c(3, 4, 12)               # frequencies Option B
#' n <- c(20, 20, 20)             # number of choices
#' strategy_postprob(k, n, list(strat1, strat2, baseline))
#' @export
strategy_postprob <- function (k, n, strategies, cores = 1){

  if (!is.null(dim(k))){
    if (is.null(dim(n)) || is.table(n))
      n <- matrix(n, nrow(k), length(n), byrow = TRUE)
    if (is.table(k))
      k <- unclass(k)
    if (cores > 1){
      cl <- makeCluster(cores)
      pp <- clusterMap(cl, strategy_postprob,
                       k = as.list(data.frame(t(k))),
                       n = as.list(data.frame(t(n))),
                       MoreArgs = list(strategies = strategies),
                       SIMPLIFY = TRUE, .scheduling = "dynamic")
      stopCluster(cl)
    } else {
      pp <- mapply(strategy_postprob,
                   k = as.list(data.frame(t(k))),
                   n = as.list(data.frame(t(n))),
                   MoreArgs = list(strategies = strategies))
    }
    if (is.matrix(pp)){
      pp <- t(pp)
      rownames(pp) <- rownames(k)
    }

  } else {
    marginal <- sapply(strategies, function(strat)
      strategy_marginal(k, n, strat))
    pp <- model_weights(marginal)
  }
  pp
}
