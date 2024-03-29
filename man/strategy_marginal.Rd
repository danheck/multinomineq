% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strategy_marginal.R
\name{strategy_marginal}
\alias{strategy_marginal}
\title{Log-Marginal Likelihood for Decision Strategy}
\usage{
strategy_marginal(k, n, strategy)
}
\arguments{
\item{k}{observed frequencies of Option B.
Either a vector or a matrix/data frame (one person per row).}

\item{n}{vector with the number of choices per item type.}

\item{strategy}{a list that defines the predictions of a strategy, see\code{\link{strategy_multiattribute}}.}
}
\description{
Computes the logarithm of the marginal likelihood, defined as the integral
over the likelihood function weighted by the prior distribution of the error probabilities.
}
\examples{
k <- c(1, 11, 18)
n <- c(20, 20, 20)
# pattern: A, A, B with constant error e<.50
strat <- list(
  pattern = c(-1, -1, 1),
  c = .5, ordered = FALSE,
  prior = c(1, 1)
)
m1 <- strategy_marginal(k, n, strat)
m1

# pattern: A, B, B with ordered error e1<e3<e2<.50
strat2 <- list(
  pattern = c(-1, 3, 2),
  c = .5, ordered = TRUE,
  prior = c(1, 1)
)
m2 <- strategy_marginal(k, n, strat2)
m2

# Bayes factor: Model 2 vs. Model 1
exp(m2 - m1)
}
