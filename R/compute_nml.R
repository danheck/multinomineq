#########################
# Compute NML complexity for decision strategies
#########################

# -log(a(theta)); luckiness not normalized
# => complexity of order-constrained models can be compared
# default:    standard NML; luck=c(1,1)
# uniform BF: inverse of Jeffreys prior; luck=c(1.5, 1.5)
luckiness <- function(par, luck = c(1, 1)){
  dd <- 0
  if (length(par) > 0)
    dd <- sum(dbeta(par, luck[1], luck[2], log = TRUE))
  dd
}

ll_luckiness <- function(par, k, n, prediction, luck = c(1,1)){
  loglik(par, k, n, prediction) + luckiness(par, luck)
}

#' Maximum-likelihood Estimate
#'
#' Get ML estimate (possible weighted by luckinesss function).
#'
#' @inheritParams compute_cnml
#' @param k vector with observed frequencies for Option B
#' @template details_prediction
#' @export
maximize_ll <- function(k, n, prediction, c = .5,
                        luck = c(1,1), n.fit = 5){

  # analytic ML estimate, boundary correction and order constraints
  est <- estimate_par(k, n, prediction, luck = luck)
  start <- adjust_par(est, prediction, c)
  npar <- length(est)

  # check for unordered models (EQW, TTB, WADD, baseline)
  # oo <- optim(ll_luckiness, lower=0, upper=c, control = list(fn=-1),
  #                k = k, n = n, prediction=prediction, luck = luck)
  if (npar > 1 && !is.null(attr(prediction, "ordered"))){
    start <- adjust_par(est, prediction, c, BOUND)
    ui <- rbind(diag(npar), 0) - rbind(0, diag(npar))
    ci <- c(rep(0, npar), -c)
    tryCatch(oo <- constrOptim(start, ll_luckiness, grad = NULL,
                               k = k, n = n, prediction = prediction, luck = luck,
                               ui = ui, ci = ci, control = list(fnscale = -1)),
             error = function(e) {print(e);
               cat("\n\n  (optimization failed with start values =", start, ")\n")})
    cnt <- 1
    while (cnt < n.fit){
      start <- sort(runif(npar, 0, c))
      oo2 <- constrOptim(start, ll_luckiness, grad = NULL,
                         k = k, n = n, prediction = prediction, luck = luck,
                         ui = ui, ci = ci, control = list(fnscale = -1))
      oo2
      cnt <- cnt + 1
      if(oo2$value > oo$value) oo <- oo2
    }
    start <- oo$par
  }
  ll <- ll_luckiness(start, k, n, prediction, luck)
  list(loglik = ll, est = start, luck = luck)
}

#' Compute NML Complexity
#'
#' Enumerates discrete data space and sums the maximum likelihood values of all data vectors.
#'
#' @param prediction a vector of strategy predictions. See details and \code{\link{predict_multiattribute}}
#' @param n vector with the number of choices per item type
#' @param c upper threshold of probabilities
#' @param luck parameters of luckiness function of LNML (i.e., parameters of a beta distribution). \code{luck = c(1.5, 1.5)} is equivalent to uniform prior for the Bayes factor
#' @param n.fit number of repeated fitting runs with random starting values
#' @param cores number of processing units to be used
#' @template details_prediction
#' @examples
#' # predict: A/A/B with error probabilities e2,e3,e1<.20
#' p <- c(-3, -1, +2)
#' compute_cnml(p, n = c(2, 2, 2), c = .20)
#'
#' # predict: A/A/B with error probabilities e2<e3<e1<.50
#' p <- c(-3, -1, +2)
#' attr(p, "ordered") <- TRUE
#' compute_cnml(p, n = c(2, 2, 2), c = .50)
#' @export
compute_cnml <- function(prediction, n, c = .50, luck = c(1, 1),
                         n.fit = 3, cores = 1){
  if (is.list(prediction)){
    if (length(c) == 1)
      c <- rep(c, length(prediction))
    res <- mapply(compute_cnml, prediction = prediction, c = c,
                  MoreArgs = list(n=n, luck=luck, n.fit=n.fit, cores=cores),
                  SIMPLIFY = FALSE)
  } else {
    check_knp(n, n, prediction)
    check_luck(luck)

    # generate all possible data and fit models
    dat <- do.call("expand.grid", lapply(n, function(nn) 0:nn))
    time <- system.time({
      if (cores > 1){
        cl <- makeCluster(cores)
        clusterExport(cl, c( "luck", "n.fit", "c", "n"),
                      envir=environment())
        lls <- parSapplyLB(cl, as.list(data.frame(t(dat))),
                           # lls <- parApply(cl = cl, X = dat, MARGIN = 1,
                           function(xx) maximize_ll(k = xx, n = n, prediction = prediction,
                                                    luck = luck, c = c,
                                                    n.fit = n.fit)$loglik)
        stopCluster(cl)
      } else {
        lls <- apply(X = dat, MARGIN = 1,
                     function(xx) maximize_ll(k = xx, n = n, prediction = prediction,
                                              luck = luck, c = c,
                                              n.fit = n.fit)$loglik)
      }
    })

    cnml <- log(sum(exp(lls)))
    res <- list(cnml = cnml, n = n, prediction = prediction, c = c,
                luck = luck, time = time["elapsed"])
  }
  res
}


#' Strategy Selection Using NML
#'
#' Compute (L)NML for observed frequencies.
#' @inheritParams compute_cnml
#' @param k observed frequencies of Option B.
#'          Either a vector or a matrix/data frame (one person per row)
#' @param cnml NML complexity values computed with \code{\link{compute_cnml}}.
#'             Note that the sample size \code{n} must match.
#' @export
select_nml <- function(k, n, cnml, n.fit = 5, cores = 1){

  if (is.matrix(k) || is.data.frame(k)){
    if (cores > 1){
      cl <- makeCluster(cores)
      nml <- parApply(cl = cl, k, 1, select_nml, n = n, cnml = cnml, n.fit = n.fit)
      stopCluster(cl)
    } else {
      nml <- apply(k, 1, select_nml, n = n, cnml = cnml, n.fit = n.fit)
    }
    if (is.matrix(nml)) nml <- t(nml)

  } else {
  k <- unlist(k)
  n <- unlist(n)
  check_knp(k, n, n)
  check_cnml(cnml, n)
  lls <- sapply(cnml, function(cc)
    maximize_ll(k = k, n = n, prediction = cc$prediction, c = cc$c,
                luck = cc$luck, n.fit = n.fit)$loglik)

  nml <- - lls + sapply(cnml, "[[", "cnml")
  strat_names <- sapply(cnml, function(cc) attr(cc, "label"))
  sel <- !sapply(strat_names, is.null)
  names(nml)[sel] <- strat_names[sel]
  }
  nml
}
