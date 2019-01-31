#include <RcppArmadillo.h>
#include <functions.h>
#include <progress.hpp>

using namespace Rcpp;

// sample from truncated gamma distribution using the inverse cdf method
// [[Rcpp::export]]
double rgamma_trunc(const double shape, const double rate,
                    const double min, const double max){
  if (min >= max)
    stop("Error in truncated gamma: Truncation boundaries violate min<max!\n");
  double pmin = R::pgamma(min, shape, 1/rate, 0, false);  // scale = 1/rate
  double pmax = R::pgamma(max, shape, 1/rate, 0, false);
  double u = R::runif(0, 1);
  return R::qgamma(pmin + u * (pmax - pmin), shape, 1/rate, 0, false);
}

// sample from truncated beta distribution using the inverse cdf method
// [[Rcpp::export]]
double rbeta_trunc(const double shape1, const double shape2,
                   const double min, const double max){
  if (min >= max)
    stop("Error in truncated beta: Truncation boundaries violate min<max!\n");
  double pmin = R::pbeta(min, shape1, shape2, 0, false);
  double pmax = R::pbeta(max, shape1, shape2, 0, false);
  double u = R::runif(0, 1);
  return R::qbeta(pmin + u * (pmax - pmin), shape1, shape2, 0, false);
}

// beta-distribution sampling (for conjugate beta)
arma::mat rbeta_mat(const unsigned int n, const arma::vec shape1, const arma::vec shape2){
  const unsigned int D = shape1.n_elem;
  mat X(n, D);
  for (unsigned int d = 0 ; d < D ; d++){
    X.col(d) = vec(rbeta(n, shape1(d), shape2(d)));
  }
  return X;
}

arma::mat rbeta_mat(const unsigned int n, const unsigned int D,
                    const double shape1, const double shape2){
  return rbeta_mat(n, shape1 * ones(D), shape2 * ones(D));
}

arma::ivec rpb_vec(const arma::vec prob, const arma::vec n){
  unsigned int I = prob.n_elem;
  arma::ivec k(I);
  for(unsigned int i = 0; i < I; i++)
    k(i) = R::rbinom(n(i), prob(i));
  return k;
}

// [[Rcpp::export]]
NumericVector ppp_bin(const arma::mat& prob, const arma::vec& k, const arma::vec& n){
  const unsigned int M = prob.n_rows;
  vec tt, kpp, x2o(M), x2p(M);
  for(unsigned int m = 0; m < M; m++){
    tt = conv_to< colvec >::from(prob.row(m));
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
arma::mat sampling_bin(const arma::vec& k, const arma::vec& n,
                       const arma::mat& A, const arma::vec& b,
                       const arma::vec& prior, const unsigned int M, arma::vec start,
                       const unsigned int burnin = 5, const bool progress = true){
  unsigned int D = A.n_cols;    // dimensions
  mat X(D, M + burnin);     // initialize posterior (column-major ordering)
  X.col(0) = start_random(A, b, M, start);

  Progress p(M, progress);
  double bmax, bmin;
  vec rhs = b;
  uvec Apos, Aneg;
  IntegerVector idx = seq_len(D) - 1;
  unsigned int j = 0;
  for (unsigned int i = 1 ; i < M + burnin; i++){
    p.increment();   // update progress bar
    if(i % 1000 == 0) Rcpp::checkUserInterrupt();
    // copy old and update to new parameters:
    X.col(i) = X.col(i-1);
    idx = sample(idx, D, false);

    for (unsigned int m = 0; m < D; m++){
      j = idx(m);
      // get min/max for truncated beta:
      bmax = 1.;
      bmin = 0.;
      rhs = (b - A * X.col(i) + A.col(j) * X(j,i))/ A.col(j);
      Aneg = find(A.col(j) < 0);
      Apos = find(A.col(j) > 0);
      if (!Aneg.is_empty())
        bmin = fmax(0., rhs(Aneg).max());
      if (!Apos.is_empty())
        bmax = fmin(1, rhs(Apos).min());
      X(j,i) = rbeta_trunc(k(j) + prior(0), n(j) - k(j) + prior(1), bmin, bmax);
    }
  }
  X.shed_cols(0, burnin - 1);
  return X.t();
}


// [[Rcpp::export]]
NumericMatrix count_bin(const arma::vec& k, const arma::vec& n,
                        const arma::mat& A, const arma::vec& b,
                        const arma::vec& prior, const unsigned int M,
                        const unsigned int batch, const bool progress = true){
  Progress p(M/batch, progress);
  int count = 0, todo = M;
  mat X;
  while (todo > 0){
    p.increment();   // update progress bar
    Rcpp::checkUserInterrupt();
    // count prior and posterior samples that match constraints:
    X = rbeta_mat(fmin(todo,batch), k + prior(0), n - k + prior(1));
    count = count  + count_samples(X, A, b);
    todo = todo - batch;
  }
  return results(count, M, A.n_rows);
}



// count samples stepwise to get volume of polytope
// [[Rcpp::export]]
NumericMatrix count_stepwise_bin(const arma::vec& k, arma::vec& n,
                                 const arma::mat& A, arma::vec& b, const arma::vec& prior,
                                 arma::vec M, arma::vec steps, const unsigned int batch,
                                 arma::vec start, const unsigned int burnin,
                                 const bool progress = true){
  steps = sort_steps(steps - 1, A.n_rows);  // C++ --> R indexing!!
  unsigned int S = steps.n_elem;    // number of unique steps
  if (M.n_elem == 1)
    M = M(0) * ones(S);

  // dynamic start values for each step
  mat starts(steps.n_elem, A.n_cols);
  for (unsigned int s = 0; s < steps.n_elem; s++)
    starts.row(s) = start.t();

  mat sample;
  uvec inside_idx;
  vec inside, count = zeros(S);
  count(0) = count_bin(k, n, A.rows(0, steps(0)),
        b.subvec(0, steps(0)), prior, M(0), batch, false)(0,0);

  // go from  A[0:steps(s-1),] ---> A[0:steps(s),]
  for (unsigned int s = 1; s < S ; s++){
    Rcpp::checkUserInterrupt();
    if (progress) Rcout << (s==1 ? " step: " : " , ") << s;
    sample = sampling_bin(k, n, A.rows(0, steps(s-1)), b.subvec(0, steps(s-1)),
                          prior, M(s), starts.row(s-1).t(), burnin, false);
    inside = inside_Ab(sample, A.rows(steps(s-1) + 1, steps(s)),
                       b.subvec(steps(s-1) + 1, steps(s)));
    if (any(inside)){
      inside_idx = find(inside);
      starts.row(s) = sample.row(inside_idx(inside_idx.n_elem - 1));
    }
    count(s) = accu(inside);
  }
  if (progress) Rcout << "\n";
  return results(count, M, steps + 1); // C++ --> R indexing
}

// [[Rcpp::export]]
NumericMatrix count_auto_bin(const arma::vec& k, const arma::vec& n,
                             const arma::mat& A, const arma::vec& b, const arma::vec& prior,
                             arma::vec count, arma::vec M, arma::vec steps,
                             const unsigned int M_iter, const unsigned int cmin,
                             const unsigned int maxiter, arma::vec start,
                             const unsigned int burnin, const bool progress = true){
  steps = steps - 1; // R --> C++ indexing
  vec inside;
  uvec inside_idx;
  mat prob;

  // dynamic starting values for each step
  mat starts(steps.n_elem, A.n_cols);
  for (unsigned int s = 0; s < steps.n_elem; s++)
    starts.row(s) = start.t();

  int i, from, iter = 0;  // from: can be negative, thus not unsigned!
  while (count.min() < cmin){
    Rcpp::checkUserInterrupt();
    if (progress && iter % (maxiter/100) == 0)
      Rcout << (iter== 0 ? " current cmin: " : " , ") << count.min();
    iter += 1;
    // find step with smallest count:
    i = count.index_min();
    from = (i == 0) ? -1 : steps(i-1);  // i==0: start at A(0,:)

    // sampling for step with minimal counts
    if (i == 0){
      prob = rbeta_mat(M_iter, k + prior(0), n - k + prior(1));
    }  else {
      prob = sampling_bin(k, n, A.rows(0, from), b.subvec(0, from),
                           prior, M_iter, starts.row(i - 1).t(), burnin, false);
      starts.row(i - 1) = prob.row(M_iter - 1);
    }
    // count samples
    inside = inside_Ab(prob, A.rows(from + 1, steps(i)),
                       b.subvec(from + 1, steps(i)));
    if (any(inside)){
      inside_idx = find(inside);
      starts.row(i) = prob.row(inside_idx(inside_idx.n_elem - 1));
    }
    count(i) += accu(inside);
    M(i) += M_iter;
    // TODO: compute precision intead of cmin
  }
  if (progress) Rcout << "\n";
  return results(count, M, steps + 1); // C++ --> R indexing
}


// ------ hit-and-run: does not seem to increase efficiency
// not as efficient as random-direction Gibbs sampling
//  (requires constraints 0<p<1  in Ab representation)
// [[Rcpp::export]]
arma::mat sampling_hitandrun(const arma::mat& A, const arma::vec& b,
                             const unsigned int M, arma::vec start,
                             const unsigned int burnin = 5, const bool progress = true){
  unsigned int D = A.n_cols;    // dimensions
  mat X(D, M + burnin);     // initialize posterior (column-major ordering)
  X.col(0) = start_random(A, b, M, start);

  Progress p(M, progress);
  double bmax = 1, bmin = 0;
  vec rhs, u, x, z;
  uvec Apos, Aneg;
  for (unsigned int i = 1 ; i < M + burnin; i++){
    p.increment();   // update progress bar
    if(i % 1000 == 0) Rcpp::checkUserInterrupt();
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
  X.shed_cols(0, burnin - 1);
  return X.t();
}
