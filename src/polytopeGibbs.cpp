# include <RcppArmadillo.h>
// # include <omp.h>
// / [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

// sample from truncated beta distribution using the inverse cdf method
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
  arma::mat X(n, D);
  for (int d = 0 ; d < D ; d++)
  {
    X.col(d) = arma::vec(rbeta(n, shape1(d), shape2(d)));
  }
  return X;
}

arma::mat rbeta_mat(int n, int D, double shape1, double shape2)
{
  return rbeta_mat(n, shape1 * arma::ones(D), shape2 * arma::ones(D));
}

// count number of samples that adhere to constraint A*x <= b
// X: samples (rows: replications; cols: D dimensions)
// [[Rcpp::export]]
int count_samples(arma::mat X, arma::mat A, arma::vec b)
{
  int cnt = 0, idx = 0;
  bool inside;
  for (int m = 0 ; m < X.n_rows ; m++)
  {
    inside = true;
    idx = 0;
    while (inside && idx < b.n_elem)
    {
      inside = arma::dot(X.row(m), A.row(idx)) <= b(idx);
      idx += 1;
    }
    cnt = cnt + (int)inside;
  }
  return cnt;
}

// DIFFICULT:
// arma::vec get_interior_point(arma::mat A, arma::vec b)
// {
//   arma::vec u = arma::randu(A.n_rows);
//
//   return x;
// }

// posterior sampling for polytope: A*x <= b
// => uniform prior sampling if  k=n=b(0,...,0)
// start: permissible starting values (randomly drawn if start[1]==-1)
// [[Rcpp::export]]
arma::mat sampling_binary_cpp(arma::vec k, arma::vec n,
                              arma::mat A, arma::vec b, arma::vec prior,
                              int M, arma::vec start)
{
  const int K = A.n_rows, D = A.n_cols;  // vertices + dimensions
  arma::mat spost(M, D);     // initialize posterior

  // find permissible starting values:
  int cnt = 0;
  if (start(0) == -1)
  {
    bool search = true;
    while (search && cnt < fmax(M, 1000))
    {
      cnt++;
      start.randu(D);
      search = any(A * start > b);
    }
    if (cnt == fmax(M, 1000))
      stop("Could not find starting values within the polytope");
  }
  spost.row(0) = start.t();

  // initialize boundaries for truncated beta sampling:
  double bnd, bmax, bmin;
  for (unsigned int i = 1 ; i < M ; i++)
  {
    // copy old values:
    spost.row(i) = spost.row(i-1);
    // update parameters:
    for (int j = 0 ; j < A.n_cols ; j++)
    {
      // initialize default boundaries:
      bmax = 1.; bmin = 0.;
      // get min/max for truncated beta:   [TODO: vectorize!]
      for (int v = 0 ; v < K ; v++)
      {
        bnd = b(v) - arma::dot(A.row(v),  spost.row(i)) + A(v,j) * spost(i,j);
        if (A(v,j) < 0)
          bmin = fmax( bnd/A(v,j), bmin);
        else if (A(v,j) > 0)
          bmax = fmin(bnd/A(v,j), bmax);
      }
      spost(i,j) =
        rbeta_trunc(k(j) + prior(0), n(j) - k(j) + prior(1), bmin, bmax);
    }
  }
  return spost;
}


// [[Rcpp::export]]
NumericVector count_polytope_cpp(arma::vec k, arma::vec n,
                                 arma::mat A, arma::vec b, arma::vec prior,
                                 int M, int batch = 10000)
{
  int count = 0, todo = M;
  arma::mat X;
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
  steps = arma::resize(steps, steps.n_elem + 1, 1);
  steps(steps.n_elem - 1) = max;
  return arma::sort(arma::unique(steps - 1));  // C++ indexing
}

// count samples to get volume of polytope
// (splits M samples into batches of size "batch" to decrease memory usage)
// [[Rcpp::export]]
List count_stepwise(arma::vec k, arma::vec n,
                    arma::mat A, arma::vec b, arma::vec prior,
                    arma::vec M, arma::vec steps, int batch = 5000)
{
  steps = sort_steps(steps, A.n_rows);
  int S = steps.n_elem;    // number of unique steps
  int D = A.n_cols;        // number of dimensions/parameters
  if (M.n_elem == 1)
    M = M(0) * arma::ones(S);
  arma::mat sample;
  arma::vec cnt = arma::zeros(S);

  int todo = M(0);
  while (todo > 0)
  {
    // sampling from unit cube and go to first polytope
    sample = rbeta_mat(fmin(todo,batch), k + prior(0), n - k + prior(1));
    cnt(0) = cnt(0) +
      count_samples(sample, A.rows(0, steps(0)), b.subvec(0, steps(0)));
    todo = todo - batch;
  }

  for (int s = 1; s < S ; s++)
  {
    sample =
      sampling_binary_cpp(k, n,
                          A.rows(0, steps(s-1)),
                          b.subvec(0, steps(s-1)),
                          prior, (int)M(s), -arma::ones(1));
    cnt(s) = cnt(s) +
      count_samples(sample,
                    A.rows(steps(s-1) + 1, steps(s)),
                    b.subvec(steps(s-1) + 1, steps(s)));
  }
  return Rcpp::List::create(Named("integral") = prod(cnt/M),
                            Named("count") = as<NumericVector>(wrap(cnt)),
                            Named("M") = as<NumericVector>(wrap(M)),
                            Named("steps") = steps + 1);  // R indexing
}

