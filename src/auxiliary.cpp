#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector adj_iterative(NumericVector par,
                            const double c = .50,
                            const double DIFF_BOUND = 0.0) {

  int I = par.length();
  par[I-1] = par[I-1] > c - DIFF_BOUND ? c - DIFF_BOUND : par[I-1];
  for (int i = I-2; i >= 0; i--)
  {
    par[i] = par[i] <= par[i+1] ? par[i] : par[i+1];
  }
  return par;
}

