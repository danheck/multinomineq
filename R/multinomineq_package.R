#' multinomineq: Bayesian Inference for Inequality-Constrained Multinomial Models
#'
#' Implements Gibbs sampling and Bayes factors for multinomial models with
#' convex, linear-inequality constraints on the probability parameters. This
#' includes models that predict a linear order of binomial probabilities
#' (e.g., p1 < p2 < p3 < .50) and mixture models, which assume that the
#' parameter vector p must be inside the convex hull of a finite number of
#' vertices.
#'
#' Inequality-constrained multinomial models have applications in the
#' area of judgment and decision making to fit and test random utility models
#' (Regenwetter, M., Dana, J., & Davis-Stober, C. P. (2011).
#'   Transitivity of preferences. Psychological Review, 118, 42–56,
#'   <doi:10.1037/a0021150>) or to perform outcome-based strategy classification
#' (i.e., to select the strategy that provides the best account for a vector of
#'   observed choice frequencies; Heck, D. W., Hilbig, B. E., & Moshagen, M.
#'   (2017). From information processing to decisions: Formalizing and comparing
#'   probabilistic choice models. Cognitive Psychology, 96, 26–40.
#'   <doi:10.1016/j.cogpsych.2017.05.003>).
#'
#' The method is implemented for multiple purposes:
#' \itemize{
#'  \item Risky decisions between different gambles to test choice axioms
#'        such as transitivity (Regenwetter et al., 2012, 2014).
#'        See \code{\link{regenwetter2012}}
#'  \item Testing deterministic axioms of measurement and choice (Karabatsos, 2005; Myung et al., 2005).
#'        See \code{\link{bf_nonlinear}}
#'  \item Multiattribute decisions for probabilistic inferences involving strategies such as
#'        Take-the-best (TTB) vs. weighted additive (WADD; Bröder & Schiffer, 2003; Heck et al., 2017)
#'        See \code{\link{heck2017}} and \code{\link{hilbig2014}}
#'  \item Fitting and testing nonparametric item response theory models (Karabatsos & Sheu, 2004).
#'        See \code{\link{karabatsos2004}}
#'  \item Statistical inference for order-constrained contingency tables (Klugkist et al., 2007, 2010).
#'        See \code{\link{bf_nonlinear}}
#'  \item Testing stochastic dominance of response time distributions (Heathcote et al., 2010).
#'        See \code{\link{stochdom_bf}}
#' }
#'
#' For convex polytopes, the transformation of vertex to inequality representation
#' requires the R package \code{rPorta} available at \url{https://github.com/TasCL/rPorta}
#'
#' @author Daniel W. Heck
#' @docType package
#'
#' @importFrom stats sd dbeta rbinom dbinom pbeta integrate runif constrOptim rbeta rmultinom
#' @importFrom parallel parSapply clusterExport makeCluster stopCluster parApply parSapplyLB clusterMap
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom utils head combn
#' @importFrom graphics hist
#' @importFrom quadprog solve.QP
#' @useDynLib "multinomineq", .registration=TRUE
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
#' @template ref_hoijtink2014
#' @template ref_klugkist2007
#' @template ref_klugkist2010
#' @template ref_heathcote2010
"_PACKAGE"


