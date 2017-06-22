#########################
#' Compute Log-Marginal Probabilities for Strategies
#'
#' Computes the logarithm of the integral over the likelihood function weighted by the prior distribution of the parameters.
#'
#' @param prior prior parameters of beta distribution, theta[i]~dbeta(prior[1], prior[2])
#' @inheritParams compute_cnml
#' @inheritParams select_nml
#' @export
compute_marginal <- function(b, n, prediction, prior = c(1, 1), c=0.5){

  n_par <- get_par_number(prediction)
  adherence <- estimate_par(b, n, - prediction, prob = FALSE)
  error <- estimate_par(b, n, prediction, prob = FALSE)
  pm <- 0

  if (c == 1 && n_par == length(b) && is.null(attr(prediction, "ordered"))){
    # baseline
    pm <- sum(lbeta(b + prior[1], n - b + prior[2]) - lbeta(prior[1], prior[2]))

  } else if (n_par == 1){
    # deterministic (constant error parameter)
    s1 <- sum(n)-adherence + prior[1]
    s2 <- adherence + prior[2]
    pm <- pbeta(c, s1, s2, log = TRUE) + lbeta(s1, s2) + # prior*lik (unnormalized)
      - pbeta(c, prior[1], prior[2], log = TRUE) - lbeta(prior[1], prior[2])
    # (normalization of order-constrained prior)

  } else if (n_par > 1){
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
    pm <- log(integrate(ff, 0, c, idx = n_par)$value) +
      lfactorial(n_par) - n_par* (pbeta(c, prior[1], prior[2], log = TRUE) +
                                    lbeta(prior[1], prior[2]))
  }

  if (any(prediction == 0)){
    # GUESS
    pm <- pm + sum(n[prediction == 0]) * log(.50)
  }
  sum(lchoose(n, b)) + pm
}


#' Strategy Selection Using the Bayes Factor
#'
#' @inheritParams compute_cnml
#' @inheritParams select_nml
#' @inheritParams compute_marginal
#' @export
select_bf <- function(b, n, prediction, c = .5, prior = c(1, 1), ...){
  UseMethod("select_bf", b)
}

#' @rdname select_bf
#' @export
select_bf.default <- function (b, n, prediction, c = .5, prior = c(1, 1)){

  if (is.numeric(prediction)){
    marginal <- compute_marginal(b, n, prediction, prior, c)
  } else {
    marginal <- sapply(prediction, function(p)
      compute_marginal(b, n, p, prior, c))
  }
  marginal
}

#' @rdname select_bf
#' @export
select_bf.matrix <- function (b, n, prediction, c = .5,
                              prior = c(1, 1), cores = 1){
  if (cores > 1){
    cl <- makeCluster(cores)
    marg <- parApply(cl, b, 1, select_bf, n = n, prediction = prediction,
                     c = c, prior = prior)
    stopCluster(cl)
  } else {
    marg <- apply(b, 1, select_bf, n = n, prediction = prediction,
                  c = c, prior = prior)
  }
  t(marg)
}

#' @rdname select_bf
#' @export
select_bf.data.frame <- function (b, n, prediction, c = .5,
                                  prior = c(1, 1), cores = 1){
  b <- as.matrix(b)
  UseMethod("select_bf", b)
}
