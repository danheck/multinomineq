#' multinomineq: Bayesian Inference for Inequality-Constrained Multinomial Models
#'
#' @description{
#' \if{html}{\figure{multinomineq.png}{options: width=120 alt ="multinomineq" style='float: right'}}
#' \if{latex}{\figure{multinomineq.png}{options: width=0.5in}}
#'
#' Implements Gibbs sampling and Bayes factors for multinomial models with
#' convex, linear-inequality constraints on the probability parameters. This
#' includes models that predict a linear order of binomial probabilities
#' (e.g., p1 < p2 < p3 < .50) and mixture models, which assume that the
#' parameter vector p must be inside the convex hull of a finite number of
#' vertices.
#' }
#'
#' @details
#' A formal definition of inequality-constrained multinomial models and the
#' implemented computational methods for Bayesian inference is provided in:
#' \itemize{
#' \item Heck, D. W., & Davis-Stober, C. P. (2019).
#'   Multinomial models with linear inequality constraints:
#'   Overview and improvements of computational methods for Bayesian inference.
#'   Manuscript under revision. \url{https://arxiv.org/abs/1808.07140}
#' }
#'
#' Inequality-constrained multinomial models have applications in multiple areas
#' in psychology, judgement and decision making, and beyond:
#' \itemize{
#'  \item Testing choice axioms such as transitivity and random utility theory
#'        (Regenwetter et al., 2012, 2014). See \code{\link{regenwetter2012}}
#'  \item Testing deterministic axioms of measurement and choice (Karabatsos, 2005; Myung et al., 2005).
#'  \item Multiattribute decisions for probabilistic inferences involving strategies such as
#'        Take-the-best (TTB) vs. weighted additive (WADD; Bröder & Schiffer, 2003; Heck et al., 2017)
#'        See \code{\link{heck2017}} and \code{\link{hilbig2014}}
#'  \item Fitting and testing nonparametric item response theory models (Karabatsos & Sheu, 2004).
#'        See \code{\link{karabatsos2004}}
#'  \item Statistical inference for order-constrained contingency tables (Klugkist et al., 2007, 2010).
#'        See \code{\link{bf_nonlinear}}
#'  \item Testing stochastic dominance of response time distributions (Heathcote et al., 2010).
#'        See \code{\link{stochdom_bf}}
#'  \item Cognitive diagnostic assessment (Hoijtink et al., 2014).
#' }
#'
##### For convex polytopes, the transformation of the vertex to the inequality representation
##### was possible with the R package \code{rPorta} (see \url{https://github.com/TasCL/rPorta}).
##### However, \code{rPorta} cannot be compiled anymore with R>=4.0.0.
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


