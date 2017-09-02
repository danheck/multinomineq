#' stratsel: Bayesian Inference for Inequality/Order-Constrained Multinomial Models
#'
#' Strategy classification is used to select the decision strategy that
#' provides the best account for a vector of observed choice frequencies.
#' By expressing choice strategies as statistical models for binary choices,
#' model-selection methods such as the Bayes factors and minimum description length
#' can be used.
#'
#' The method is implemented for multiple purposes:
#' \itemize{
#' \item Multiattribute decisions for probabilistic inferences involving strategies such as
#'       Take-the-best (TTB) vs. weighted additive (WADD; Bröder, 2003; Heck et al., 2017)
#' \item Risky decisions between different gambles to test choice axioms
#'       such as transitivity (Regenwetter et al., 2012, 2014).
#'  \item Testing deterministic axioms of measurement and choice (Karabatsos, 2005; Myung et al., 2005).
#'  \item Fitting and testing nonparametric item response theory models (Karabatsos, 2002).
#'  \item Statistical inference for order-constrained contingency tables (Klugkist et al., 2007, 2010).
#' }
#'
#' Detailed explanations and examples can be found in the package vignette, accessible via \code{vignette("stratsel")}.
#' For convex polytopes, the transformation of vertex to inequality representation requires the R package \code{rPorta} available at \url{https://github.com/TasCL/rPorta}
#'
#' @author Daniel W. Heck
#' @docType package
#' @importFrom stats sd dbeta rbinom dbinom pbeta integrate runif constrOptim rbeta rmultinom
#' @importFrom parallel parSapply clusterExport makeCluster stopCluster parApply parSapplyLB clusterMap
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom utils head combn
#' @importFrom graphics hist
#' @importFrom quadprog solve.QP
#' @useDynLib "stratsel", .registration=TRUE
#'
#' @template ref_broder2003
#' @template ref_broder2003jepg
#' @template ref_heck2017
#' @template ref_hilbig2014
#' @template ref_regenwetter2012
#' @template ref_regenwetter2014
#' @template ref_karabatsos2005
#' @template ref_myung2005
#' @template ref_karabatsos2004
#' @template ref_hoijtink2011
#' @template ref_klugkist2007
#' @template ref_klugkist2010

"_PACKAGE"


