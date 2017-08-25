#include <RcppArmadillo.h>
#include <functions.h>
#include <progress.hpp>

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

arma::ivec rpb_vec(arma::vec theta, arma::vec n)
{
  int I = theta.n_elem;
  arma::ivec k(I);
  for(int i = 0; i < I; i++)
    k(i) = R::rbinom(n(i), theta(i));
  return k;
}

// [[Rcpp::export]]
NumericVector ppp_bin(arma::mat theta, arma::vec k, arma::vec n)
{
  int M = theta.n_rows;
  vec tt, kpp, x2o(M), x2p(M);
  for(int m = 0; m < M; m++)
  {
    tt = conv_to< colvec >::from(theta.row(m));
    kpp = conv_to< vec >::from(rpb_vec(tt, n));
    x2o(m) = x2(k, tt % n);
    x2p(m) = x2(kpp, tt % n);
  }
  double ppp = double(accu(x2o <= x2p)) / M;
  return NumericVector::create(Named("X2_obs") = mean(x2o),
                               Named("X2_pred") = mean(x2p),
                               Named("ppp") = ppp);
}


// posterior sampling for polytope: A*x <= b
// => uniform prior sampling if  k=n=b(0,...,0)
// start: permissible starting values (randomly drawn if start[1]==-1)
// [[Rcpp::export]]
arma::mat sampling_binomial_cpp(arma::vec k, arma::vec n,
                                arma::mat A, arma::vec b,
                                arma::vec prior, int M, arma::vec start,
                                int burnin = 5, double progress = true)
{
  int D = A.n_cols;    // dimensions
  mat X(D, M + burnin);     // initialize posterior (column-major ordering)
  X.col(0) = start_random(A, b, M, start);

  Progress p(M, progress);
  bool run = true;
  double bmax, bmin;
  vec rhs = b;
  uvec Apos, Aneg;
  IntegerVector idx = seq_len(D) - 1;
  int j=0, steps=100;
  for (int i = 1 ; i < M + burnin; i++)
  {
    p.increment();   // update progress bar
    if(run && i % steps == 0) run = !Progress::check_abort();
    if (run)
    {
      // copy old and update to new parameters:
      X.col(i) = X.col(i-1);
      idx = sample(idx, D, false);

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
  }
  X.shed_cols(0, burnin - 1);
  return X.t();
}


// [[Rcpp::export]]
NumericVector count_binomial_cpp(arma::vec k, arma::vec n,
                                 arma::mat A, arma::vec b,
                                 arma::vec prior, int M,
                                 int batch, bool progress = true)
{
  Progress p(M/batch, progress);
  bool run = true;
  int count = 0, todo = M;
  mat X;
  while (todo > 0)
  {
    p.increment();   // update progress bar
    if (run)
    {
      run = !Progress::check_abort();
      // count prior and posterior samples that match constraints:
      X = rbeta_mat(fmin(todo,batch), k + prior(0), n - k + prior(1));
      count = count  + count_samples(X, A, b);
      todo = todo - batch;
    }
  }
  return NumericVector::create(Named("integral") = (double)count / M,
                               // Named("results") = results(count, M));
                               Named("count") = count,
                               Named("M") = M);
}

// add the last order constraint and sort "steps" vector
arma::vec sort_steps(arma::vec steps, int total)
{
  steps = resize(steps, steps.n_elem + 1, 1);
  steps(steps.n_elem - 1) = total - 1;
  return sort(unique(steps));
}

// go from  A[0:from,] ---> A[0:to,]  // C++ indexing!
// [[Rcpp::export]]
int count_step(arma::vec k, arma::vec n,
               arma::mat A, arma::vec b, arma::vec prior,
               int M, int from, int to, arma::vec start, bool progress = true)
{
  mat sample =
    sampling_binomial_cpp(k, n, A.rows(0, from), b.subvec(0, from),
                          prior, M, start, 10, progress);
  return count_samples(sample, A.rows(from + 1, to), b.subvec(from + 1, to));
}

// count samples stepwise to get volume of polytope
// [[Rcpp::export]]
List count_stepwise(arma::vec k, arma::vec n,
                    arma::mat A, arma::vec b, arma::vec prior,
                    arma::vec M, arma::vec steps, int batch,
                    arma::vec start, bool progress = true)
{
  steps = sort_steps(steps, A.n_rows);  // C++ indexing!!
  int S = steps.n_elem;    // number of unique steps
  int D = A.n_cols;        // number of dimensions/parameters
  if (M.n_elem == 1)
    M = M(0) * ones(S);
  mat sample;

  vec cnt = zeros(S);
  cnt(0) = count_binomial_cpp(k, n, A.rows(0, steps(0)),
      b.subvec(0, steps(0)), prior, M(0), batch, progress)["count"];

  for (int s = 1; s < S ; s++)
    cnt(s) = count_step(k, n, A, b, prior, M(s), steps(s-1), steps(s), start, progress);

  return List::create(Named("integral") = prod(cnt / M.rows(0, S-1)),
                      // Named("results") = results(steps, cnt, M),
                      Named("count") = as<NumericVector>(wrap(cnt)),
                      Named("M") = as<NumericVector>(wrap(M)),
                      Named("steps") = steps);  // R indexing
}


// not as efficient as random-direction Gibbs sampling
//  (requires constraints 0<p<1  in Ab representation)
// [[Rcpp::export]]
arma::mat sampling_hitandrun(arma::mat A, arma::vec b, int M, arma::vec start,
                             int burnin = 5, double progress = true)
{
  int D = A.n_cols;    // dimensions
  mat X(D, M + burnin);     // initialize posterior (column-major ordering)
  X.col(0) = start_random(A, b, M, start);

  Progress p(M, progress);
  bool run = true;
  double bmax = 1, bmin = 0;
  vec rhs, u, x, z;
  uvec Apos, Aneg;
  int steps=100;
  for (int i = 1 ; i < M + burnin; i++)
  {
    p.increment();   // update progress bar
    if(run && i % steps == 0) run = !Progress::check_abort();
    if (run)
    {
      // random direction
      u = normalise(randn(D,1));
      z = A * u;
      x = X.col(i-1);
      rhs = (b - A * x) / z;

      // get min/max for truncated uniform:
      bmax = 1.; bmin = 0.;
      if (any(z < 0))
        bmin = rhs(find(z < 0)).max();
      if (any(z > 0))
        bmax = rhs(find(z > 0)).min();
      X.col(i) = x + R::runif(bmin, bmax) * u;
    }
  }
  X.shed_cols(0, burnin - 1);
  return X.t();
}
