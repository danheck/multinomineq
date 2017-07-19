#include <RcppArmadillo.h>
#include <functions.h>
using namespace Rcpp;


// sample from truncated beta distribution using the inverse cdf method
// [[Rcpp::export]]
double rbeta_trunc(double shape1, double shape2,
                   double min, double max)
{
  double pmin = R::pbeta(min, shape1, shape2, 0, false);
  double pmax = R::pbeta(max, shape1, shape2, 0, false);
  double u = R::runif(0, 1);
  return R::qbeta(pmin + u * (pmax - pmin), shape1, shape2, 0, false);
}

// beta-distribution sampling (for conjugate beta)
arma::mat rbeta_mat(int n, arma::vec shape1, arma::vec shape2)
{
  int D = shape1.n_elem;
  mat X(n, D);
  for (int d = 0 ; d < D ; d++)
  {
    X.col(d) = vec(rbeta(n, shape1(d), shape2(d)));
  }
  return X;
}

arma::mat rbeta_mat(int n, int D, double shape1, double shape2)
{
  return rbeta_mat(n, shape1 * ones(D), shape2 * ones(D));
}

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
  // int cnt = 0, idx = 0;
  // bool inside;
  // for (int m = 0 ; m < X.n_rows ; m++)
  // {
  //   inside = true;
  //   idx = 0;
  //   while (inside && idx < b.n_elem)
  //   {
  //     inside = dot(X.row(m), A.row(idx)) <= b(idx);
  //     idx += 1;
  //   }
  //   cnt = cnt + (int)inside;
  // }
  return accu(inside_Ab(X, A, b));
}

int count_samples(arma::vec X, arma::mat A, arma::vec b)
{
  mat Xm = reshape(X, 1, X.n_elem);
  return count_samples(Xm, A, b);
}


// DIFFICULT:
// arma::vec get_interior_point(arma::mat A, arma::vec b)
// {
//   arma::vec u = arma::randu(A.n_rows);
//
//   return x;
// }
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

// posterior sampling for polytope: A*x <= b
// => uniform prior sampling if  k=n=b(0,...,0)
// start: permissible starting values (randomly drawn if start[1]==-1)
// [[Rcpp::export]]
arma::mat sampling_binomial_cpp(arma::vec k, arma::vec n,
                                arma::mat A, arma::vec b,
                                arma::vec prior, int M, arma::vec start,
                                int burnin = 5)
{
  int D = A.n_cols;    // dimensions
  mat X(D, M + burnin);     // initialize posterior (column-major ordering)
  X.col(0) = start_random(A, b, M, start);
  // X.col(0) = start;

  double bmax, bmin;
  vec rhs = b;
  uvec Apos, Aneg;
  ivec idx;
  int j;
  for (int i = 1 ; i < M + burnin; i++)
  {
    // copy old and update to new parameters:
    X.col(i) = X.col(i-1);
    idx = randi(D, distr_param(0,D - 1));
    for (int m = 0; m < D; m++)
    {
      j = idx(m);
      // get min/max for truncated beta:
      bmax = 1.; bmin = 0.;
      rhs = (b - A * X.col(i) + A.col(j) * X(j,i))/ A.col(j);
      Aneg = find(A.col(j) < 0);
      Apos = find(A.col(j) > 0);
      if (!Aneg.is_empty())
        bmin = rhs(Aneg).max();
      if (!Apos.is_empty())
        bmax = rhs(Apos).min();
      X(j,i) = rbeta_trunc(k(j) + prior(0), n(j) - k(j) + prior(1), bmin, bmax);
    }
  }
  X.shed_cols(0, burnin - 1);
  return X.t();
}


// [[Rcpp::export]]
NumericVector count_binomial_cpp(arma::vec k, arma::vec n,
                                 arma::mat A, arma::vec b,
                                 arma::vec prior, int M, int batch)
{
  int count = 0, todo = M;
  mat X;
  while (todo > 0)
  {
    // count prior and posterior samples that match constraints:
    X = rbeta_mat(fmin(todo,batch), k + prior(0), n - k + prior(1));
    count = count  + count_samples(X, A, b);
    todo = todo - batch;
  }
  return NumericVector::create(Named("integral") = (double)count / M,
                               Named("count") = count,
                               Named("M") = M);
}

// add the last order constraint and sort "steps" vector
arma::vec sort_steps(arma::vec steps, int max)
{
  steps = resize(steps, steps.n_elem + 1, 1);
  steps(steps.n_elem - 1) = max;
  return sort(unique(steps - 1));  // C++ indexing
}

// count samples to get volume of polytope
// (splits M samples into batches of size "batch" to decrease memory usage)
// [[Rcpp::export]]
List count_stepwise(arma::vec k, arma::vec n,
                    arma::mat A, arma::vec b, arma::vec prior,
                    arma::vec M, arma::vec steps, int batch,
                    arma::vec start)
{
  steps = sort_steps(steps, A.n_rows);
  int S = steps.n_elem;    // number of unique steps
  int D = A.n_cols;        // number of dimensions/parameters
  if (M.n_elem == 1)
    M = M(0) * ones(S);
  mat sample;

  vec cnt = zeros(S);
  cnt(0) = count_binomial_cpp(k, n,
      A.rows(0, steps(0)),
      b.subvec(0, steps(0)), prior, M(0), batch)["count"];

  for (int s = 1; s < S ; s++)
  {
    sample =
      sampling_binomial_cpp(k, n, A.rows(0, steps(s-1)), b.subvec(0, steps(s-1)),
                            prior, (int)M(s), start);
    cnt(s) = cnt(s) +
      count_samples(sample,
                    A.rows(steps(s-1) + 1, steps(s)),
                    b.subvec(steps(s-1) + 1, steps(s)));
  }
  return List::create(Named("integral") = prod(cnt/M),
                      Named("count") = as<NumericVector>(wrap(cnt)),
                      Named("M") = as<NumericVector>(wrap(M)),
                      Named("steps") = steps + 1);  // R indexing
}

