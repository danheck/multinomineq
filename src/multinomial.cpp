#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/rmultinom.h>
#include <functions.h>
#include <progress.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat rdirichlet(int n, arma::vec alpha)
{
  int I = alpha.n_elem;
  mat X(n, I);
  double ai;
  for (int i = 0; i < I; i++)
  {
    ai = as_scalar(alpha(i));
    X.col(i) = vec(rgamma(n, ai, 1));
  }
  X = X / (sum(X, 1) * ones(1, I));
  return X;
}

//' Random Samples from the Product-Dirichlet Distribution
//'
//' Random samples from the prior/posterior (i.e., product-Dirichlet) of the unconstrained
//' product-multinomial model (the encompassing model).
//'
//' @param n number of samples
//' @param alpha Dirichlet parameters concatenated across independent conditions
//'     (e.g., a1,a2,a3,  b1,b2,b3, ..)
//' @param options the number of choice options per item type, e.g., \code{c(3,2)}
//'     for a ternary and binary condition. The sum of \code{options} must be equal to the length of \code{alpha}.
//' @examples
//' # standard uniform Dirichlet
//' rpdirichlet(5, c(1,1,1,1), 4)
//'
//' # two ternary outcomes: (a1,a2,a3,  b1,b2,b3)
//' rpdirichlet(5, c(9,5,1,  3,6,6), c(3,3))
//' @export
// [[Rcpp::export]]
arma::mat rpdirichlet(int n, arma::vec alpha, arma::vec options)
{
  int D = options.n_elem;  // number of conditions/item types
  mat X(n, alpha.n_elem);
  vec sel = cumsum(options) - 1;
  sel.insert_rows(0, - ones(1, 1));
  int k,l;
  for (int i = 0; i < D; i++)
  {
    k = sel(i) + 1;
    l = sel(i + 1);
    X.cols(k, l) = rdirichlet(n, alpha.rows(k, l));
  }
  return X;
}

// omit last column for each item type
//  [a1 a2 a3  b1 b2]  => [a1 a2  b1]
// [[Rcpp::export]]
arma::mat rpdirichlet_free(int n, arma::vec alpha, arma::vec options)
{
  mat X = rpdirichlet(n, alpha, options);
  int D = options.n_elem;
  vec dep_idx = cumsum(options) - 1;
  for (int i = D - 1; i >= 0; i--)
  {
    X.shed_col(dep_idx(i));
  }
  return X;
}

// remove dependent categories:  (a1,a2, b1,b2,b3)  =>  (a1,  b1,b2)
// [[Rcpp::export]]
arma::vec shed_options(arma::vec x, arma::vec options)
{
  int D = options.n_elem;
  vec sel = cumsum(options) - 1;
  for (int i = D - 1; i >= 0; i--)
    x.shed_row(sel(i));
  return x;
}

// replicate option per vector:  c(10,20)  =>  c(10,10,  20,20,20)
// [[Rcpp::export]]
arma::vec rep_options(arma::vec x, arma::vec options)
{
  vec y;
  for (int i = 0; i < options.n_elem; i++)
    y = join_cols(y,repmat(x(i) * ones(1,1), options(i), 1));
  return y;
}

// sum per option, length equal to options (for k or prior)
arma::vec sum_options_short(arma::vec k, arma::vec options)
{
  vec n = zeros(options.n_elem, 1);
  int o = 0, cnt = 0;
  for (int i = 0; i < k.n_elem; i++)
  {
    if (cnt < options(o))
    {
      cnt += 1; // still within the same option
    } else {
      o += 1;   // new option
      cnt = 0;
    }
    n(o) += k(i);
  }
  return n;
}

// sum per option, length equal to input (for k or prior)
// [[Rcpp::export]]
arma::vec sum_options(arma::vec k, arma::vec options)
{
  vec n = sum_options_short(k, options);
  return rep_options(n, options);
}

arma::ivec rpm_vec(arma::vec theta, arma::vec n, arma::vec options)
{
  int I = options.n_elem, s0 = 0, s1, nn;
  ivec k, tmp;
  vec tt;
  for(int i = 0; i < I; i++)
  {
    s1 = s0 + options(i) - 1;
    tt = theta.rows(s0, s1);
    s0 = s1 + 1;
    nn = as_scalar(n(i));
    tmp = Rcpp::RcppArmadillo::rmultinom(nn, as<NumericVector>(wrap(tt)));
    k = join_cols(k, tmp);
  }
  return k;
}


// [[Rcpp::export]]
arma::imat rpm_mat(arma::mat theta, arma::vec n, arma::vec options)
{
  imat k(theta.n_cols, theta.n_rows);
  for(int i = 0; i < theta.n_rows; i++)
    k.col(i) = rpm_vec(theta.row(i).t(), n, options);
  return k.t();
}

// theta: with fixed probabiltiies!
// [[Rcpp::export]]
NumericVector ppp_mult(arma::mat theta, arma::vec k, arma::vec options)
{
  int M = theta.n_rows;
  vec n = sum_options(k, options);
  vec n_short = sum_options_short(k, options);
  vec tt, kpp, x2o(M), x2p(M);
  for(int m = 0; m < M; m++)
  {
    tt = conv_to< colvec >::from(theta.row(m));
    kpp = conv_to< vec >::from(rpm_vec(tt, n_short, options));
    x2o(m) = x2(k, tt % n);
    x2p(m) = x2(kpp, tt % n);
  }
  double ppp = double(accu(x2o <= x2p)) / M;
  return NumericVector::create(Named("X2_obs") = mean(x2o),
                               Named("X2_pred") = mean(x2p),
                               Named("ppp") = ppp);
}

// [[Rcpp::export]]
arma::mat sampling_mult(arma::vec k, arma::vec options,
                        arma::mat A, arma::vec b,
                        arma::vec prior, int M, arma::vec start,
                        int burnin = 5, double progress = true)
{
  int D = A.n_cols;  // dimensions
  mat X(D, M + burnin);       // initialize posterior, column-major order
  X.col(0) = start_random(A, b, M, start);

  // last options for each multinomial (not included in matrix A)
  uvec idx_J = conv_to<uvec>::from(cumsum(options) - 1);
  vec beta_J = rep_options(k(idx_J) + prior(idx_J), options - 1);
  vec beta_j = shed_options(k + prior, options);
  // lower/upper index within each multinomial:
  vec j_up  = rep_options(cumsum(options - 1) - 1, options - 1);
  vec j_low = j_up - rep_options(options - 2, options - 1);

  Progress p(M, progress);
  bool run = true;
  double bmax, bmin, s;
  vec rhs = b;
  uvec Apos, Aneg;
  IntegerVector idx = seq_len(D) - 1;
  int j=0, tmp=100;
  for (int i = 1 ; i < M  + burnin; i++)
  {
    p.increment();   // update progress bar
    if(run && i % tmp == 0) run = !Progress::check_abort();
    if (run)
    {
      X.col(i) = X.col(i-1);
      idx = sample(idx, D, false);
      for (int m = 0; m < D; m++)
      {
        j = m; //idx(m);
        // scaled truncated beta:
        s = 1 - accu(X.col(i).rows(j_low(j),j_up(j))) + X(j,i);  // scaling within multinomial
        bmax = 1; bmin = 0;
        rhs = (b - A * X.col(i) + A.col(j) * X(j,i)) / (A.col(j) * s);
        Aneg = find(A.col(j) < 0);
        Apos = find(A.col(j) > 0);
        if (!Aneg.is_empty())
          bmin = rhs(Aneg).max();
        if (!Apos.is_empty())
          bmax = rhs(Apos).min();
        X(j,i) = s * rbeta_trunc(beta_j(j), beta_J(j), bmin, bmax);
      }
    }
  }
  X.shed_cols(0,burnin - 1);
  return X.t();
}


// [[Rcpp::export]]
NumericMatrix count_mult(arma::vec k, arma::vec options,
                         arma::mat A, arma::vec b,
                         arma::vec prior, int M,
                         int batch, bool progress = true)
{
  Progress p(M/batch, progress);
  bool run = true;
  int count = 0, todo = M;
  mat X(batch, k.n_elem);
  while (todo > 0)
  {
    p.increment();   // update progress bar
    if (run)
    {
      // count prior and posterior samples that match constraints:
      X = rpdirichlet_free(fmin(todo,batch), k + prior, options);
      count = count  + count_samples(X, A, b);
      todo = todo - batch;
    }
  }
  return results(count, M, A.n_rows);
}

// go from  A[0:from,] ---> A[0:to,]  // C++ indexing!
// [[Rcpp::export]]
int count_step_mult(arma::vec k, arma::vec options,
                    arma::mat A, arma::vec b, arma::vec prior,
                    int M, int from, int to, arma::vec start, bool progress = true)
{
  mat sample = sampling_mult(k, options, A.rows(0, from), b.subvec(0, from),
                             prior, M, start, 10, progress);
  return count_samples(sample, A.rows(from + 1, to), b.subvec(from + 1, to));
}

// [[Rcpp::export]]
NumericMatrix count_stepwise_multi(arma::vec k, arma::vec options,
                                   arma::mat A, arma::vec b, arma::vec prior,
                                   arma::vec M, arma::vec steps, int batch,
                                   arma::vec start, bool progress = true)
{
  steps = sort_steps(steps - 1, A.n_rows);  // R --> C++ indexing
  int S = steps.n_elem;    // number of unique steps
  int D = A.n_cols;        // number of dimensions/parameters
  if (M.n_elem == 1)
    M = M(0) * ones(S);

  vec count = zeros(S);
  count(0) = count_mult(k, options, A.rows(0, steps(0)),
        b.subvec(0, steps(0)), prior, (int)M(0), batch, false)(0,0);

  for (int s = 1; s < S ; s++)
  {
    Rcpp::checkUserInterrupt();
    if (progress) Rcout << (s==1 ? " step: " : " , ") << s;
    count(s) = count_step_mult(k, options, A, b, prior, M(s),
          steps(s-1), steps(s), start, false);
  }
  if (progress) Rcout << "\n";
  return results(count, M, steps + 1); // C++ --> R indexing
}

// [[Rcpp::export]]
NumericMatrix count_auto_mult(arma::vec k, arma::vec options,
                              arma::mat A, arma::vec b, arma::vec prior,
                              arma::vec count, arma::vec M, arma::vec steps,
                              int M_iter, int cmin, int maxiter,
                              arma::vec start, bool progress = true)
{
  steps = steps - 1; // R --> C++ indexing
  vec inside;
  mat theta;
  mat starts(steps.n_elem, A.n_cols);
  for (int s = 0; s < steps.n_elem; s++) starts.row(s) = start.t(); // dynamic start values
  int i, from, iter=0;
  while(count.min() < cmin)
  {
    Rcpp::checkUserInterrupt();
    if (progress && iter % (maxiter/100) == 0)
      Rcout << (iter== 0 ? " current cmin: " : " , ") << count.min();
    iter += 1;
    // find step with smallest count:
    i = count.index_min();
    from = (i == 0) ? -1 : steps(i-1);  // i==0: start at A(0,:)

    // sampling for step with minimal counts
    if (i == 0)
    {
      theta = rpdirichlet_free(M_iter, k + prior, options);
    }  else {
      theta = sampling_mult(k, options, A.rows(0, from), b.subvec(0, from),
                            prior, M_iter, starts.row(i).t(), 5, false);
      starts.row(i) = theta.row(M_iter - 1);
    }
    // count samples
    inside = inside_Ab(theta, A.rows(from + 1, steps(i)),
                       b.subvec(from + 1, steps(i)));
    count(i) += accu(inside);
    M(i) += M_iter;
    // store last successful sample (update starting value)
    // compute precision
    // Rcout << count << "\n##\n";
  }
  if (progress) Rcout << "\n";
  return results(count, M, steps + 1); // C++ --> R indexing
}
