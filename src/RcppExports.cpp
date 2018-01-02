// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// inside_Ab
arma::vec inside_Ab(const arma::mat& X, const arma::mat& A, const arma::vec& b);
RcppExport SEXP _multinomineq_inside_Ab(SEXP XSEXP, SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(inside_Ab(X, A, b));
    return rcpp_result_gen;
END_RCPP
}
// count_samples
int count_samples(const arma::mat& X, const arma::mat& A, const arma::vec& b);
RcppExport SEXP _multinomineq_count_samples(SEXP XSEXP, SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(count_samples(X, A, b));
    return rcpp_result_gen;
END_RCPP
}
// start_random
arma::vec start_random(const arma::mat& A, const arma::vec& b, const unsigned int M, arma::vec start);
RcppExport SEXP _multinomineq_start_random(SEXP ASEXP, SEXP bSEXP, SEXP MSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(start_random(A, b, M, start));
    return rcpp_result_gen;
END_RCPP
}
// adj_iterative
NumericVector adj_iterative(NumericVector par, const double c, const double DIFF_BOUND);
RcppExport SEXP _multinomineq_adj_iterative(SEXP parSEXP, SEXP cSEXP, SEXP DIFF_BOUNDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type DIFF_BOUND(DIFF_BOUNDSEXP);
    rcpp_result_gen = Rcpp::wrap(adj_iterative(par, c, DIFF_BOUND));
    return rcpp_result_gen;
END_RCPP
}
// rbeta_trunc
double rbeta_trunc(const double shape1, const double shape2, const double min, const double max);
RcppExport SEXP _multinomineq_rbeta_trunc(SEXP shape1SEXP, SEXP shape2SEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< const double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< const double >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(rbeta_trunc(shape1, shape2, min, max));
    return rcpp_result_gen;
END_RCPP
}
// ppp_bin
NumericVector ppp_bin(const arma::mat& prob, const arma::vec& k, const arma::vec& n);
RcppExport SEXP _multinomineq_ppp_bin(SEXP probSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(ppp_bin(prob, k, n));
    return rcpp_result_gen;
END_RCPP
}
// sampling_bin
arma::mat sampling_bin(const arma::vec& k, const arma::vec& n, const arma::mat& A, const arma::vec& b, const arma::vec& prior, const unsigned int M, arma::vec start, const unsigned int burnin, const bool progress);
RcppExport SEXP _multinomineq_sampling_bin(SEXP kSEXP, SEXP nSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP startSEXP, SEXP burninSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_bin(k, n, A, b, prior, M, start, burnin, progress));
    return rcpp_result_gen;
END_RCPP
}
// count_bin
NumericMatrix count_bin(const arma::vec& k, const arma::vec& n, const arma::mat& A, const arma::vec& b, const arma::vec& prior, const unsigned int M, const unsigned int batch, const bool progress);
RcppExport SEXP _multinomineq_count_bin(SEXP kSEXP, SEXP nSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP batchSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type batch(batchSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(count_bin(k, n, A, b, prior, M, batch, progress));
    return rcpp_result_gen;
END_RCPP
}
// count_stepwise_bin
NumericMatrix count_stepwise_bin(const arma::vec& k, arma::vec& n, const arma::mat& A, arma::vec& b, const arma::vec& prior, arma::vec M, arma::vec steps, const unsigned int batch, arma::vec start, const unsigned int burnin, const bool progress);
RcppExport SEXP _multinomineq_count_stepwise_bin(SEXP kSEXP, SEXP nSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP stepsSEXP, SEXP batchSEXP, SEXP startSEXP, SEXP burninSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type batch(batchSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(count_stepwise_bin(k, n, A, b, prior, M, steps, batch, start, burnin, progress));
    return rcpp_result_gen;
END_RCPP
}
// count_auto_bin
NumericMatrix count_auto_bin(const arma::vec& k, const arma::vec& n, const arma::mat& A, const arma::vec& b, const arma::vec& prior, arma::vec count, arma::vec M, arma::vec steps, const unsigned int M_iter, const unsigned int cmin, const unsigned int maxiter, arma::vec start, const unsigned int burnin, const bool progress);
RcppExport SEXP _multinomineq_count_auto_bin(SEXP kSEXP, SEXP nSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP countSEXP, SEXP MSEXP, SEXP stepsSEXP, SEXP M_iterSEXP, SEXP cminSEXP, SEXP maxiterSEXP, SEXP startSEXP, SEXP burninSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type count(countSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type M_iter(M_iterSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type cmin(cminSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(count_auto_bin(k, n, A, b, prior, count, M, steps, M_iter, cmin, maxiter, start, burnin, progress));
    return rcpp_result_gen;
END_RCPP
}
// sampling_hitandrun
arma::mat sampling_hitandrun(const arma::mat& A, const arma::vec& b, const unsigned int M, arma::vec start, const unsigned int burnin, const bool progress);
RcppExport SEXP _multinomineq_sampling_hitandrun(SEXP ASEXP, SEXP bSEXP, SEXP MSEXP, SEXP startSEXP, SEXP burninSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_hitandrun(A, b, M, start, burnin, progress));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet
arma::mat rdirichlet(const unsigned int n, const arma::vec alpha);
RcppExport SEXP _multinomineq_rdirichlet(SEXP nSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet(n, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rpdirichlet
arma::mat rpdirichlet(const unsigned int n, const arma::vec alpha, const arma::vec options, const bool p_drop);
RcppExport SEXP _multinomineq_rpdirichlet(SEXP nSEXP, SEXP alphaSEXP, SEXP optionsSEXP, SEXP p_dropSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< const bool >::type p_drop(p_dropSEXP);
    rcpp_result_gen = Rcpp::wrap(rpdirichlet(n, alpha, options, p_drop));
    return rcpp_result_gen;
END_RCPP
}
// rep_options
arma::vec rep_options(arma::vec x, const arma::vec options);
RcppExport SEXP _multinomineq_rep_options(SEXP xSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(rep_options(x, options));
    return rcpp_result_gen;
END_RCPP
}
// sum_options
arma::vec sum_options(const arma::vec k, const arma::vec options);
RcppExport SEXP _multinomineq_sum_options(SEXP kSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_options(k, options));
    return rcpp_result_gen;
END_RCPP
}
// rpm_mat
arma::imat rpm_mat(const arma::mat& prob, const arma::vec& n, const arma::vec& options);
RcppExport SEXP _multinomineq_rpm_mat(SEXP probSEXP, SEXP nSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(rpm_mat(prob, n, options));
    return rcpp_result_gen;
END_RCPP
}
// ppp_mult
NumericVector ppp_mult(const arma::mat& prob, const arma::vec& k, const arma::vec& options);
RcppExport SEXP _multinomineq_ppp_mult(SEXP probSEXP, SEXP kSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(ppp_mult(prob, k, options));
    return rcpp_result_gen;
END_RCPP
}
// sampling_mult
arma::mat sampling_mult(const arma::vec& k, const arma::vec& options, const arma::mat& A, const arma::vec& b, const arma::vec& prior, const unsigned int M, arma::vec start, const unsigned int burnin, const bool progress);
RcppExport SEXP _multinomineq_sampling_mult(SEXP kSEXP, SEXP optionsSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP startSEXP, SEXP burninSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_mult(k, options, A, b, prior, M, start, burnin, progress));
    return rcpp_result_gen;
END_RCPP
}
// count_mult
NumericMatrix count_mult(const arma::vec& k, const arma::vec& options, const arma::mat& A, const arma::vec& b, const arma::vec& prior, const unsigned int M, const unsigned int batch, const bool progress);
RcppExport SEXP _multinomineq_count_mult(SEXP kSEXP, SEXP optionsSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP batchSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type batch(batchSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(count_mult(k, options, A, b, prior, M, batch, progress));
    return rcpp_result_gen;
END_RCPP
}
// count_stepwise_multi
NumericMatrix count_stepwise_multi(const arma::vec& k, const arma::vec& options, const arma::mat& A, const arma::vec& b, const arma::vec& prior, arma::vec M, arma::vec steps, const unsigned int batch, arma::vec start, const unsigned int burnin, const bool progress);
RcppExport SEXP _multinomineq_count_stepwise_multi(SEXP kSEXP, SEXP optionsSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP stepsSEXP, SEXP batchSEXP, SEXP startSEXP, SEXP burninSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type batch(batchSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(count_stepwise_multi(k, options, A, b, prior, M, steps, batch, start, burnin, progress));
    return rcpp_result_gen;
END_RCPP
}
// count_auto_mult
NumericMatrix count_auto_mult(const arma::vec& k, const arma::vec& options, const arma::mat& A, const arma::vec& b, const arma::vec& prior, arma::vec count, arma::vec M, arma::vec steps, const unsigned int M_iter, const unsigned int cmin, const unsigned int maxiter, arma::vec start, const unsigned int burnin, const bool progress);
RcppExport SEXP _multinomineq_count_auto_mult(SEXP kSEXP, SEXP optionsSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP countSEXP, SEXP MSEXP, SEXP stepsSEXP, SEXP M_iterSEXP, SEXP cminSEXP, SEXP maxiterSEXP, SEXP startSEXP, SEXP burninSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type count(countSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type M_iter(M_iterSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type cmin(cminSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(count_auto_mult(k, options, A, b, prior, count, M, steps, M_iter, cmin, maxiter, start, burnin, progress));
    return rcpp_result_gen;
END_RCPP
}
