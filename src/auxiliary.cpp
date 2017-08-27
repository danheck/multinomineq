#include <RcppArmadillo.h>
#include <functions.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;



// count number of samples that adhere to constraint A*x <= b
// X: samples (rows: replications; cols: D dimensions)
// [[Rcpp::export]]
arma::vec inside_Ab(arma::mat X, arma::mat A, arma::vec b)
{
  vec i(X.n_rows);
  int idx = 0;
  bool inside;
  for (int m = 0 ; m < X.n_rows ; m++)
  {
    inside = true;
    idx = 0;
    while (inside && idx < b.n_elem)
    {
      inside = dot(X.row(m), A.row(idx)) <= b(idx);
      idx += 1;
    }
    i(m) = inside;
  }
  return i;
}

// [[Rcpp::export]]
int count_samples(arma::mat X, arma::mat A, arma::vec b)
{
  return accu(inside_Ab(X, A, b));
}

int count_samples(arma::vec X, arma::mat A, arma::vec b)
{
  mat Xm = reshape(X, 1, X.n_elem);
  return count_samples(Xm, A, b);
}

// find permissible starting values:
// [[Rcpp::export]]
arma::vec start_random(arma::mat A, arma::vec b, int M, arma::vec start)
{
  if (start(0) == -1)
  {
    int cnt = 0;
    bool search = true;
    while (search && cnt < fmax(M, 1000))
    {
      cnt++;
      start.randu(A.n_cols);
      search = 1 - count_samples(start, A, b);
    }
    if (cnt == fmax(M, 1000))
      stop("Could not find starting values within the polytope");
  }
  return start;
}

double x2(arma::vec o, arma::vec e){
  return accu(pow(o - e, 2) / e);
}


// get named NumericMatrix with results
NumericMatrix results(arma::vec count, arma::vec M, arma::vec steps)
{
  int S = count.n_elem;
  mat res = join_rows(count, join_rows(M.rows(0, S-1), steps));
  NumericMatrix results = as<NumericMatrix>(wrap(res));
  colnames(results) = CharacterVector::create("count", "M", "steps");
  return results;
}

NumericMatrix results(int count, int M, int steps)
{
  NumericMatrix results(1,3);
  results(0,0) = count;
  results(0,1) = M;
  results(0,2) = steps;
  colnames(results) = CharacterVector::create("count", "M", "steps");
  return results;
}

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
