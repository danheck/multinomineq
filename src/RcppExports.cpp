// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// adj_iterative
NumericVector adj_iterative(NumericVector par, const double c, const double DIFF_BOUND);
RcppExport SEXP stratsel_adj_iterative(SEXP parSEXP, SEXP cSEXP, SEXP DIFF_BOUNDSEXP) {
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
double rbeta_trunc(double shape1, double shape2, double min, double max);
RcppExport SEXP stratsel_rbeta_trunc(SEXP shape1SEXP, SEXP shape2SEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< double >::type min(minSEXP);
    Rcpp::traits::input_parameter< double >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(rbeta_trunc(shape1, shape2, min, max));
    return rcpp_result_gen;
END_RCPP
}
// count_samples
int count_samples(arma::mat X, arma::mat A, arma::vec b);
RcppExport SEXP stratsel_count_samples(SEXP XSEXP, SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(count_samples(X, A, b));
    return rcpp_result_gen;
END_RCPP
}
// start_random
arma::vec start_random(arma::mat A, arma::vec b, int M, arma::vec start);
RcppExport SEXP stratsel_start_random(SEXP ASEXP, SEXP bSEXP, SEXP MSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(start_random(A, b, M, start));
    return rcpp_result_gen;
END_RCPP
}
// sampling_binomial_cpp
arma::mat sampling_binomial_cpp(arma::mat A, arma::vec b, arma::vec k, arma::vec n, arma::vec prior, int M, arma::vec start, int burnin);
RcppExport SEXP stratsel_sampling_binomial_cpp(SEXP ASEXP, SEXP bSEXP, SEXP kSEXP, SEXP nSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP startSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_binomial_cpp(A, b, k, n, prior, M, start, burnin));
    return rcpp_result_gen;
END_RCPP
}
// count_binomial_cpp
NumericVector count_binomial_cpp(arma::mat A, arma::vec b, arma::vec k, arma::vec n, arma::vec prior, int M, int batch);
RcppExport SEXP stratsel_count_binomial_cpp(SEXP ASEXP, SEXP bSEXP, SEXP kSEXP, SEXP nSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP batchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type batch(batchSEXP);
    rcpp_result_gen = Rcpp::wrap(count_binomial_cpp(A, b, k, n, prior, M, batch));
    return rcpp_result_gen;
END_RCPP
}
// count_stepwise
List count_stepwise(arma::mat A, arma::vec b, arma::vec k, arma::vec n, arma::vec prior, arma::vec M, arma::vec steps, int batch);
RcppExport SEXP stratsel_count_stepwise(SEXP ASEXP, SEXP bSEXP, SEXP kSEXP, SEXP nSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP stepsSEXP, SEXP batchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< int >::type batch(batchSEXP);
    rcpp_result_gen = Rcpp::wrap(count_stepwise(A, b, k, n, prior, M, steps, batch));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet
arma::mat rdirichlet(int n, arma::vec alpha);
RcppExport SEXP stratsel_rdirichlet(SEXP nSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet(n, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rpdirichlet
arma::mat rpdirichlet(int n, arma::vec alpha, arma::vec options);
RcppExport SEXP stratsel_rpdirichlet(SEXP nSEXP, SEXP alphaSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(rpdirichlet(n, alpha, options));
    return rcpp_result_gen;
END_RCPP
}
// rpdirichlet_free
arma::mat rpdirichlet_free(int n, arma::vec alpha, arma::vec options);
RcppExport SEXP stratsel_rpdirichlet_free(SEXP nSEXP, SEXP alphaSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(rpdirichlet_free(n, alpha, options));
    return rcpp_result_gen;
END_RCPP
}
// count_multinomial_cpp
NumericVector count_multinomial_cpp(arma::mat A, arma::vec b, arma::vec options, arma::vec k, arma::vec prior, int M, int batch);
RcppExport SEXP stratsel_count_multinomial_cpp(SEXP ASEXP, SEXP bSEXP, SEXP optionsSEXP, SEXP kSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP batchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type batch(batchSEXP);
    rcpp_result_gen = Rcpp::wrap(count_multinomial_cpp(A, b, options, k, prior, M, batch));
    return rcpp_result_gen;
END_RCPP
}
// shed_options
arma::vec shed_options(arma::vec x, arma::vec options);
RcppExport SEXP stratsel_shed_options(SEXP xSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(shed_options(x, options));
    return rcpp_result_gen;
END_RCPP
}
// rep_options
arma::vec rep_options(arma::vec x, arma::vec options);
RcppExport SEXP stratsel_rep_options(SEXP xSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(rep_options(x, options));
    return rcpp_result_gen;
END_RCPP
}
// sum_options
arma::vec sum_options(arma::vec k, arma::vec options);
RcppExport SEXP stratsel_sum_options(SEXP kSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_options(k, options));
    return rcpp_result_gen;
END_RCPP
}
// sampling_multinomial_cpp
arma::mat sampling_multinomial_cpp(arma::mat A, arma::vec b, arma::vec options, arma::vec k, arma::vec prior, int M, arma::vec start, int burnin);
RcppExport SEXP stratsel_sampling_multinomial_cpp(SEXP ASEXP, SEXP bSEXP, SEXP optionsSEXP, SEXP kSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP startSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_multinomial_cpp(A, b, options, k, prior, M, start, burnin));
    return rcpp_result_gen;
END_RCPP
}
