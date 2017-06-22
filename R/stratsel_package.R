#' stratsel:
#'
#' Outcome-based
#'
#' Detailed explanations and examples can be found in the package vignette, accessible via \code{vignette("stratsel")}
#' @author Daniel W. Heck
#' @docType package
#' @importFrom stats sd dbeta rbinom
#' @importFrom parallel parSapply clusterExport makeCluster stopCluster parApply
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib "stratsel", .registration=TRUE
#'
#' @template ref_heck17
#' @template ref_hilbig14
#' @template ref_broder03
#'
"_PACKAGE"


