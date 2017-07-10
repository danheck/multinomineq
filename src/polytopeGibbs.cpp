# include <RcppArmadillo.h>
// # include <omp.h>
using namespace Rcpp;

// sample from truncated beta distribution
double rbeta_trunc(double shape1, double shape2,
                   double min, double max)
{
  // cumulative density for lower/upper bound:
  double pmin = R::pbeta(min, shape1, shape2, 0, false);
  double pmax = R::pbeta(max, shape1, shape2, 0, false);

  // inverse cdf method:
  double u = R::runif(0, 1);
  return R::qbeta(pmin + u * (pmax - pmin), shape1, shape2, 0, false);
}

// beta-distribution sampling (for conjugate beta)
arma::mat rbeta_mat(int n, arma::vec shape1, arma::vec shape2)
{
  int D = shape1.n_elem;
  arma::mat x = arma::mat(n, D);
  for(unsigned int d = 0 ; d < D ; d++)
  {
    x.col(d) = arma::vec(rbeta(n, shape1(d), shape2(d)));
  }
  return x;
}

arma::mat rbeta_mat(int n, int D, double shape1, double shape2)
{
  return rbeta_mat(n, shape1 * arma::ones(D), shape2 * arma::ones(D));
}

// count number of samples that adhere to constraint A*x <= b
// x: samples (rows: D dimensions, cols: M replications)
// [[Rcpp::export]]
int count_samples(arma::mat x, arma::mat A, arma::vec b)
{
  int cnt = 0.;
  arma::rowvec ct = b.t();
  arma::mat At = A.t();
  for (unsigned int m = 0 ; m < x.n_rows ; m++)
  {
    cnt = cnt + all(x.row(m) * At <= ct);
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

// standard encompassing approach
// (splits M samples into batches of size "batch" to decrease memory usage)
// [[Rcpp::export]]
NumericVector encompassing_bf(arma::vec k, arma::vec n,
                              arma::mat A, arma::vec b, arma::vec prior,
                              int M, int batch = 5000)
{
  unsigned  int D = A.n_cols;
  int npost,nprior = 0;
  arma::mat sprior,spost;
  int m = M;
  while (m > 0)
  {
    // count prior and posterior samples that match constraints:
    sprior = rbeta_mat(fmin(m, batch), D, prior(0), prior(1));
    spost  = rbeta_mat(fmin(m, batch), k + prior(0), n - k + prior(1));
    npost  = npost  + count_samples(spost,  A, b);
    nprior = nprior + count_samples(sprior, A, b);
    m = m - batch;
  }
  return NumericVector::create(Named("bf") = (double)npost / nprior,
                               Named("posterior") = npost,
                               Named("prior") = nprior,
                               Named("M") = M);
}


// posterior sampling for polytope: A*x <= b
// => uniform prior sampling if  k=n=b(0,...,0)
// start: permissible starting values (randomly drawn if start[1]==-1)
// [[Rcpp::export]]
arma::mat sampling_binary_cpp(arma::vec k, arma::vec n,
                             arma::mat A, arma::vec b, arma::vec prior,
                             int M, arma::vec start)
{
  const int K = A.n_rows ;  // vertices
  const int D = A.n_cols ;  // dimensions
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

// count samples to get volume of polytope
// (splits M samples into batches of size "batch" to decrease memory usage)
// [[Rcpp::export]]
List encompassing_stepwise(arma::vec k, arma::vec n,
                           arma::mat A, arma::vec b, arma::vec prior,
                           arma::vec M, arma::vec steps, int batch = 5000){
  // int cores = omp_get_max_threads();
  // omp_set_num_threads(cores);

  steps = arma::sort(arma::unique(steps - 1));  // C++ indexing
  unsigned int D = A.n_cols;   // number of dimensions/parameters
  int S = steps.n_elem + 1;    // number of steps
  if (M.n_elem == 1)
    M = double(M(0)) * arma::ones(S);
  int m = M(0);
  arma::mat sample;
  arma::vec cnt = arma::zeros(S);

  while (m > 0)
  {
    // sampling from unit cube and go to first polytope
    sample = rbeta_mat(fmin(m,batch), k + prior(0), n - k + prior(1));
    cnt(0) = cnt(0) +
      count_samples(sample, A.rows(0, steps(0)), b.subvec(0, steps(0)));
    m = m - batch;
  }

  for (int s = 1; s < S ; s++)
  {
    // prior/posterior sampling from constrained polytope
    sample =
      sampling_binary_cpp(k, n, A.rows(0, steps(s-1)), b.subvec(0, steps(s-1)),
                         prior, int(M(s)), -arma::ones(1));
    if (s < S - 1)
      cnt(s) = cnt(s) +
        count_samples(sample, A.rows(0, steps(s)), b.subvec(0, steps(s)));
    else
      cnt(s) = cnt(s) + count_samples(sample, A, b);
  }
  return Rcpp::List::create(Named("integral") = prod(cnt/M),
                            Named("counts") = cnt,
                            Named("M") = M,
                            Named("steps") = steps + 1);  // R indexing
}

