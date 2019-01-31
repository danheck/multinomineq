#include <RcppArmadillo.h>
#include <functions.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// count number of samples that adhere to constraint A*x <= b
bool inside_Ab(const arma::vec& X, const arma::mat& A, const arma::vec& b){
  unsigned int idx = 0;
  bool inside = true;
  while (inside && idx < b.n_elem){
    inside = dot(X, A.row(idx)) <= b(idx);
    idx += 1;
  }
  return inside;
}

// X: samples (rows: replications; cols: D dimensions)
// [[Rcpp::export]]
arma::vec inside_Ab(const arma::mat& X, const arma::mat& A, const arma::vec& b){
  vec i(X.n_rows);
  for (unsigned int m = 0 ; m < X.n_rows ; m++){
    i(m) = inside_Ab(vec(X.row(m).t()), A, b);
  }
  return i;
}

// [[Rcpp::export]]
int count_samples(const arma::mat& X, const arma::mat& A, const arma::vec& b){
  return accu(inside_Ab(X, A, b));
}

int count_samples(const arma::vec& X, const arma::mat& A, const arma::vec& b){
  mat Xmat = reshape(X, 1, X.n_elem);
  return count_samples(Xmat, A, b);
}

// find permissible starting values:
// [[Rcpp::export]]
arma::vec start_random(const arma::mat& A, const arma::vec& b,
                       const unsigned int M, arma::vec start){
  if (start(0) == -1){
    unsigned int cnt = 0;
    bool inside = false;
    while (!inside && cnt < fmax(M, 5000)){
      cnt++;
      start.randu(A.n_cols);
      inside = inside_Ab(start, A, b);
    }
    if (cnt == fmax(M, 5000))
      stop("Could not find random starting value within the polytope.");
  }
  // mat prob = reshape(start, 1, start.n_elem);
  // if (!as_scalar(inside_Ab(prob, A, b))){
  if (!inside_Ab(start, A, b)){
    Rcout << "A = \n" << A << "\nb = " << b.t() << "\nstart = " << start.t();
    stop("Starting value: Does not satisfy  A*x<b.\n");
  }
  return start;
}

// add the last order constraint and sort "steps" vector (C++ indexing!)
arma::vec sort_steps(arma::vec steps, const unsigned int A_rows){
  if (steps.max() != A_rows - 1){
    steps = resize(steps, steps.n_elem + 1, 1);
    steps(steps.n_elem - 1) = A_rows - 1;
  }
  return sort(unique(steps));
}

double x2(const arma::vec& o, const arma::vec& e){
  return accu(pow(o - e, 2) / e);
}


// get named NumericMatrix with results
NumericMatrix results(const arma::vec& count, const arma::vec& M,
                      const arma::vec& steps){
  unsigned int S = count.n_elem;
  mat res = join_rows(count, join_rows(M.rows(0, S-1), steps));
  NumericMatrix results = as<NumericMatrix>(wrap(res));
  colnames(results) = CharacterVector::create("count", "M", "steps");
  return results;
}

NumericMatrix results(const unsigned int count, const unsigned int M,
                      const unsigned int steps){
  NumericMatrix results(1,3);
  results(0,0) = count;
  results(0,1) = M;
  results(0,2) = steps;
  colnames(results) = CharacterVector::create("count", "M", "steps");
  return results;
}

// for strategy models: find appropriate value inside linear order cosntraints
// [[Rcpp::export]]
NumericVector adj_iterative(NumericVector par,
                            const double c = .50,
                            const double DIFF_BOUND = 0.0) {
  const unsigned int I = par.length();
  par[I-1] = par[I-1] > c - DIFF_BOUND ? c - DIFF_BOUND : par[I-1];
  for (int i = I-2; i >= 0; i--){
    par[i] = par[i] <= par[i+1] ? par[i] : par[i+1];
  }
  return par;
}
