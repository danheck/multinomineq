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
// count_samples
int count_samples(arma::mat x, arma::mat A, arma::vec b);
RcppExport SEXP stratsel_count_samples(SEXP xSEXP, SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(count_samples(x, A, b));
    return rcpp_result_gen;
END_RCPP
}
// bf_encompassing
NumericVector bf_encompassing(arma::vec k, arma::vec n, arma::mat A, arma::vec b, arma::vec prior, int M, int batch);
RcppExport SEXP stratsel_bf_encompassing(SEXP kSEXP, SEXP nSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP batchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type batch(batchSEXP);
    rcpp_result_gen = Rcpp::wrap(bf_encompassing(k, n, A, b, prior, M, batch));
    return rcpp_result_gen;
END_RCPP
}
// sampling_posterior
arma::mat sampling_posterior(arma::vec k, arma::vec n, arma::mat A, arma::vec b, arma::vec prior, int M, arma::vec start);
RcppExport SEXP stratsel_sampling_posterior(SEXP kSEXP, SEXP nSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_posterior(k, n, A, b, prior, M, start));
    return rcpp_result_gen;
END_RCPP
}
// encompassing_stepwise
List encompassing_stepwise(arma::vec k, arma::vec n, arma::mat A, arma::vec b, arma::vec prior, arma::vec M, arma::vec steps, int batch);
RcppExport SEXP stratsel_encompassing_stepwise(SEXP kSEXP, SEXP nSEXP, SEXP ASEXP, SEXP bSEXP, SEXP priorSEXP, SEXP MSEXP, SEXP stepsSEXP, SEXP batchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< int >::type batch(batchSEXP);
    rcpp_result_gen = Rcpp::wrap(encompassing_stepwise(k, n, A, b, prior, M, steps, batch));
    return rcpp_result_gen;
END_RCPP
}
