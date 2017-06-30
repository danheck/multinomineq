#########################
#' Compute Log-Marginal Probabilities for Strategies
#'
#' Computes the logarithm of the integral over the likelihood function weighted by the prior distribution of the parameters.
#'
#' @param prior prior parameters of beta distribution, theta[i]~dbeta(prior[1], prior[2])
#' @inheritParams compute_cnml
#' @inheritParams select_nml
#' @examples
#' k <- c(1,11,18)
#' n <- c(20, 20, 20)
#' # prediction: A, A, B with constant error e<.50
#' pred <- c(-1, -1, 1)
#' m1 <- compute_marginal(k, n, pred)
#' m1
#'
#' # prediction: A, B, B with ordered error e1<e3<e2<.50
#' pred2 <- c(-1, 3, 2)
#' attr(pred2, "ordered") <- TRUE
#' m2 <- compute_marginal(k, n, pred2)
#' m2
#'
#' # Bayes factor: Model 2 vs. Model 1
#' exp(m2 - m1)
#' @export
compute_marginal <- function (k, n, prediction, c = 0.5, prior = c(1, 1)){
  check_knpcp(k, n, prediction, c, prior)

  n_par <- get_par_number(prediction)
  adherence <- estimate_par(k, n, - prediction, prob = FALSE)
  error     <- estimate_par(k, n,   prediction, prob = FALSE)
  pm <- 0

  if (c == 1 && n_par == length(k) && is.null(attr(prediction, "ordered"))){
    # baseline
    pm <- sum(lbeta(k + prior[1], n - k + prior[2]) - lbeta(prior[1], prior[2]))

  } else if (n_par == 1){
    # deterministic (constant error parameter)
    s1 <- error + prior[1]
    s2 <- adherence + prior[2]
    pm <- pbeta(c, s1, s2, log.p = TRUE) + lbeta(s1, s2) + # prior*lik (unnormalized)
      - pbeta(c, prior[1], prior[2], log.p = TRUE) - lbeta(prior[1], prior[2])
    # (normalization of order-constrained prior)

  } else if (n_par > 1 && n_par <= 6){
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
      lfactorial(n_par) - n_par* (pbeta(c, prior[1], prior[2], log.p = TRUE) +
                                    lbeta(prior[1], prior[2]))
  } else {

  }

  if (any(prediction == 0)){
    # GUESS
    pm <- pm + sum(n[prediction == 0]) * log(.50)
  }
  sum(lchoose(n, k)) + pm
}


#' Strategy Selection Using the Bayes Factor
#'
#' Computes the posterior model probabilities for multiple strategies (with equal prior model probabilities).
#'
#' @inheritParams compute_cnml
#' @inheritParams select_nml
#' @inheritParams compute_marginal
#' @param prediction.list list of model predictions.
#' @param c upper boundary for parameters/error probabilities.
#'   Note that a vector of the same length as \code{prediction.list}
#'   can be provided to use different upper bound per model.
#' @template details_prediction
#' @seealso \code{\link{compute_marginal}} and \code{\link{model_weights}}
#' @examples
#' k <- c(3, 4, 12)               # frequencies Option B
#' n <- c(20, 20, 20)             # number of choices
#' preds <- list(
#'     strategy1 = c(-1, -1, -1), # A/A/A
#'     strategy2 = c(0, 1, -1),   # guess/B/A
#'     baseline = 1:3)
#' select_bf(k, n, preds, c= c(.5, .5, 1))
#' @export
select_bf <- function (k, n, prediction.list, c = .5,
                       prior = c(1, 1), cores = 1){

  if (!is.null(dim(k))){
    if (is.null(dim(n)))
      n <- matrix(n, nrow(k), length(n), byrow = TRUE)
    if (cores > 1){
      cl <- makeCluster(cores)
      pp <- clusterMap(cl, select_bf,
                       k = as.list(data.frame(t(k))),
                       n = as.list(data.frame(t(n))),
                       MoreArgs = list(prediction.list = prediction.list,
                                       c = c, prior = prior),
                       SIMPLIFY = TRUE, .scheduling = "dynamic")
      stopCluster(cl)
    } else {
      pp <- mapply(select_bf, k = as.list(data.frame(t(k))),
                   n = as.list(data.frame(t(n))),
                   MoreArgs = list(prediction.list = prediction.list,
                                   c = c, prior = prior))
    }
    if (is.matrix(pp)){
      pp <- t(pp)
      rownames(pp) <- rownames(k)
    }

  } else {
    if (length(c) == 1)
      c <- rep(c, length(prediction.list))

    marginal <- mapply(function(p, c) compute_marginal(k, n, p, c, prior),
                       p = prediction.list, c = as.list(c))
    pp <- model_weights(marginal)
  }
  pp
}
