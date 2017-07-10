#' stratsel: Strategy Classification for Binary Choice Models
#'
#' Strategy classification is used to select the decision strategy that
#' provides the best account for a vector of observed choice frequencies.
#' By expressing choice stratgies as statistical models for binary choices,
#' model-selection methods such as the Bayes factors and minimum description length
#' can be used.
#'
#' The method is implemented for two areas:
#' \itemize{
#' \item Multiattribute decisions for probabilistic inferences involving strategies such as
#'       Take-the-best (TTB) vs. weighted additive (WADD; Bröder, 2003)
#' \item Risky decisions between different gambles to test choice axioms
#'       such as transitivity (Regenwetter, 2014).
#' }
#'
#' Detailed explanations and examples can be found in the package vignette, accessible via \code{vignette("stratsel")}.
#' For convex polytopes, the transformation of vertex to inequality representation requires the R package \code{rPorta} available at \url{https://github.com/TasCL/rPorta}
#'
#' @author Daniel W. Heck
#' @docType package
#' @importFrom stats sd dbeta rbinom dbinom pbeta integrate runif constrOptim rbeta
#' @importFrom parallel parSapply clusterExport makeCluster stopCluster parApply parSapplyLB clusterMap
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib "stratsel", .registration=TRUE
#'
#' @template ref_broder2003
#' @template ref_broder2003jepg
#' @template ref_heck2017
#' @template ref_hilbig2014
#' @template ref_karabatsos2005
#' @template ref_regenwetter2014

"_PACKAGE"


