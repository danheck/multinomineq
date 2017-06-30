#' stratsel: Strategy Classification for Binary Choice Models
#'
#' Strategy classification uses model selection methods
#' (i.e., Bayes factors and minimum description length) to select the decision strategy that
#' provides the best account for a vector of observed choice frequencies.
#' This method can be used for multiattribute decisions involving strategies such as
#' Take-the-best (TTB) vs. weighted additive (WADD; Bröder, 2003), but also for risky decisions to
#' test the transitivity of preferences and other choice axioms (Regenwetter, 2014).
#'
#' Detailed explanations and examples can be found in the package vignette, accessible via \code{vignette("stratsel")}.
#'
#' @author Daniel W. Heck
#' @docType package
#' @importFrom stats sd dbeta rbinom dbinom pbeta integrate runif constrOptim
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


