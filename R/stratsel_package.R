#' stratsel:
#'
#' Outcome-based
#'
#' Detailed explanations and examples can be found in the package vignette, accessible via \code{vignette("stratsel")}
#' @author Daniel W. Heck
#' @docType package
#' @importFrom stats sd dbeta rbinom
#' @importFrom parallel parSapply clusterExport makeCluster stopCluster
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib "stratsel", .registration=TRUE
#'
#' @references
#' Heck, D. W., Hilbig, B. E., & Moshagen, M. (2017). From information processing to decisions: Formalizing and comparing probabilistic choice models. Cognitive Psychology, 96, 26-40. https://doi.org/10.1016/j.cogpsych.2017.05.003
#'
"_PACKAGE"


