// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

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

static const R_CallMethodDef CallEntries[] = {
    {"stratsel_adj_iterative", (DL_FUNC) &stratsel_adj_iterative, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_stratsel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}