#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/rmultinom.h>
#include <functions.h>
#include <progress.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat rdirichlet(const unsigned int n, const arma::vec alpha){
  unsigned int I = alpha.n_elem;
  mat X(n, I);
  double ai;
  for (unsigned int i = 0; i < I; i++){
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
//'     (e.g., a1,a2,  b1,b2,b3)
//' @param options the number of choice options per item type, e.g., \code{c(2,3)}
//'     for a binary and ternary condition.
//'     The sum of \code{options} must be equal to the length of \code{alpha}.
//' @param drop_fixed whether the output matrix includes the last probability for each category
//'     (which is not a free parameter since probabilities must sum to one).
//'
//' @examples
//' # standard uniform Dirichlet
//' rpdirichlet(5, c(1,1,1,1), 4)
//' rpdirichlet(5, c(1,1,1,1), 4, drop_fixed = FALSE)
//'
//' # two ternary outcomes: (a1,a2,a3,  b1,b2,b3)
//' rpdirichlet(5, c(9,5,1,  3,6,6), c(3,3))
//' rpdirichlet(5, c(9,5,1,  3,6,6), c(3,3), drop_fixed = FALSE)
//' @export
// [[Rcpp::export]]
arma::mat rpdirichlet(const unsigned int n, const arma::vec alpha,
                      const arma::vec options, const bool drop_fixed = true){
  unsigned int D = options.n_elem;  // number of conditions/item types
  mat X(n, alpha.n_elem);
  vec sel = cumsum(options);
  sel.insert_rows(0, zeros(1, 1));
  unsigned int k,l;
  for (unsigned int i = 0; i < D; i++){
    k = sel(i);
    l = sel(i + 1) - 1;
    X.cols(k, l) = rdirichlet(n, alpha.rows(k, l));
  }
  if (drop_fixed){
    // omit last column for each item type: [a1 a2 a3  b1 b2]  => [a1 a2  b1]
    vec dep_idx = cumsum(options) - 1;
    for (int i = D - 1; i >= 0; i--){
      X.shed_col(dep_idx(i));
    }
  }
  return X;
}


// remove dependent categories:  (a1,a2, b1,b2,b3)  =>  (a1,  b1,b2)
arma::vec shed_options(arma::vec x, const arma::vec options){
  unsigned int D = options.n_elem;
  vec sel = cumsum(options) - 1;
  for (int i = D - 1; i >= 0; i--)
    x.shed_row(sel(i));
  return x;
}

// replicate option per vector:  c(10,20)  =>  c(10,10,  20,20,20)
// [[Rcpp::export]]
arma::vec rep_options(arma::vec x, const arma::vec options){
  vec y;
  for (unsigned int i = 0; i < options.n_elem; i++)
    y = join_cols(y,repmat(x(i) * ones(1,1), options(i), 1));
  return y;
}

// sum per option, length equal to options (for k or prior)
arma::vec sum_options_short(const arma::vec k, const arma::vec options){
  vec n = zeros(options.n_elem, 1);
  unsigned int o = 0, cnt = 0;
  for (unsigned int i = 0; i < k.n_elem; i++){
    if (cnt < options(o)){
      cnt += 1; // still within the same option
    } else {
      o += 1;   // new option
      cnt = 1;
    }
    n(o) += k(i);
  }
  return n;
}

// sum per option, length equal to input (for k or prior)
// [[Rcpp::export]]
arma::vec sum_options(const arma::vec k, const arma::vec options){
  vec n = sum_options_short(k, options);
  return rep_options(n, options);
}

arma::ivec rpm_vec(const arma::vec& prob, const arma::vec& n,
                   const arma::vec& options){
  unsigned int I = options.n_elem, s0 = 0, s1, nn;
  ivec k, tmp;
  vec tt;
  for(unsigned int i = 0; i < I; i++){
    s1 = s0 + options(i) - 1;
    tt = prob.rows(s0, s1);
    s0 = s1 + 1;
    nn = as_scalar(n(i));
    if (nn > 0){
      tmp = Rcpp::RcppArmadillo::rmultinom(nn, as<NumericVector>(wrap(tt)));
      k = join_cols(k, tmp);
    }
  }
  return k;
}


// [[Rcpp::export]]
arma::imat rpm_mat(const arma::mat& prob, const arma::vec& n, const arma::vec& options){
  imat k(prob.n_cols, prob.n_rows);
  for(unsigned int i = 0; i < prob.n_rows; i++)
    k.col(i) = rpm_vec(prob.row(i).t(), n, options);
  return k.t();
}

// prob: with fixed probabiltiies!
// [[Rcpp::export]]
NumericVector ppp_mult(const arma::mat& prob, const arma::vec& k, const arma::vec& options){
  unsigned int M = prob.n_rows;
  vec n = sum_options(k, options);
  vec n_short = sum_options_short(k, options);
  vec tt, kpp, x2o(M), x2p(M);
  for (unsigned int m = 0; m < M; m++){
    tt = conv_to<colvec>::from(prob.row(m));
    kpp = conv_to<colvec>::from(rpm_vec(tt, n_short, options));
    x2o(m) = x2(k, tt % n);
    x2p(m) = x2(kpp, tt % n);
  }
  double ppp = double(accu(x2o <= x2p)) / M;
  return NumericVector::create(Named("X2_obs") = mean(x2o),
                               Named("X2_pred") = mean(x2p),
                               Named("ppp") = ppp);
}

// [[Rcpp::export]]
arma::mat sampling_mult(const arma::vec& k, const arma::vec& options,
                        const arma::mat& A, const arma::vec& b,
                        const arma::vec& prior, const unsigned int M, arma::vec start,
                        const unsigned int burnin = 5, const bool progress = true){
  unsigned int D = A.n_cols;  // dimensions
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
  double bmax, bmin, s;
  vec rhs = b;
  uvec Apos, Aneg;
  IntegerVector idx = seq_len(D) - 1;
  unsigned int j = 0;
  for (unsigned int i = 1 ; i < M  + burnin; i++){
    p.increment();   // update progress bar
    if(i % 50 == 0) Rcpp::checkUserInterrupt();
    X.col(i) = X.col(i-1);
    idx = sample(idx, D, false);
    for (unsigned int m = 0; m < D; m++){
      j = idx(m);
      // scaling within multinomial:
      s = 1 - accu(X.col(i).rows(j_low(j),j_up(j))) + X(j,i);
      // truncation:
      rhs = (b - A * X.col(i) + A.col(j) * X(j,i)) / A.col(j);
      Aneg = find(A.col(j) < 0);
      Apos = find(A.col(j) > 0);
      bmin = 0.;
      bmax = s;
      if (!Aneg.is_empty())
        bmin = fmax(0., rhs(Aneg).max());
      if (!Apos.is_empty())
        bmax = fmin(s, rhs(Apos).min());
      X(j,i) = s * rbeta_trunc(beta_j(j), beta_J(j), bmin/s, bmax/s);
      // //### Debugging:
      // Rcout << X.col(i).t();
      // if (as_scalar(inside_Ab(X.col(i).t(), A, b)) == 0){
      //   Rcout << "rhs(Aneg).max() =" << rhs(Aneg).max() <<"--- rhs(Apos).min()=" << rhs(Apos).min();
      //   Rcout << "==>" <<i << "//" << j << ":   s="<< s <<" ;bmin=" << bmin <<" bmax=" << bmax << "\n";
      // }
    }

  }
  X.shed_cols(0,burnin - 1);
  return X.t();
}


// [[Rcpp::export]]
NumericMatrix count_mult(const arma::vec& k, const arma::vec& options,
                         const arma::mat& A, const arma::vec& b,
                         const arma::vec& prior, const unsigned int M,
                         const unsigned int batch, const bool progress = true){
  Progress p(M/batch, progress);
  int count = 0, todo = M;
  mat X(batch, k.n_elem);
  while (todo > 0){
    p.increment();   // update progress bar
    Rcpp::checkUserInterrupt();
    // count prior and posterior samples that match constraints:
    X = rpdirichlet(fmin(todo,batch), k + prior, options, true);
    count = count  + count_samples(X, A, b);
    todo = todo - batch;
  }
  return results(count, M, A.n_rows);
}


// [[Rcpp::export]]
NumericMatrix count_stepwise_multi(const arma::vec& k, const arma::vec& options,
                                   const arma::mat& A, const arma::vec& b, const arma::vec& prior,
                                   arma::vec M, arma::vec steps, const unsigned int batch,
                                   arma::vec start, const unsigned int burnin,
                                   const bool progress = true){
  steps = sort_steps(steps - 1, A.n_rows);  // R --> C++ indexing
  unsigned int S = steps.n_elem;    // number of unique steps
  if (M.n_elem == 1)
    M = M(0) * ones(S);
  mat starts(S, A.n_cols);
  for (unsigned int s = 0; s < S; s++)
    starts.row(s) = start_random(A, b, M(0), start).t(); // different starting values for each step

  mat mcmc;
  uvec inside_idx;
  vec inside, count = zeros(S);
  count(0) = count_mult(k, options, A.rows(0, steps(0)),
        b.subvec(0, steps(0)), prior, (int)M(0), batch, false)(0,0);

  // go from  A[0:steps(s-1),] ---> A[0:steps(s),]
  for (unsigned int s = 1; s < S ; s++){
    Rcpp::checkUserInterrupt();
    if (progress) Rcout << (s==1 ? " step: " : " , ") << s;
    mcmc = sampling_mult(k, options, A.rows(0, steps(s-1)), b.subvec(0, steps(s-1)),
                         prior, M(s), starts.row(s - 1).t(), burnin, false);
    inside = inside_Ab(mcmc, A.rows(steps(s-1) + 1, steps(s)),
                       b.subvec(steps(s-1) + 1, steps(s)));
    if (any(inside)){
      inside_idx = find(inside);
      starts.row(s) = mcmc.row(inside_idx(inside_idx.n_elem - 1));
    }
    count(s) = accu(inside);
  }
  if (progress) Rcout << "\n";
  return results(count, M, steps + 1); // C++ --> R indexing
}

// [[Rcpp::export]]
NumericMatrix count_auto_mult(const arma::vec& k, const arma::vec& options,
                              const arma::mat& A, const arma::vec& b, const arma::vec& prior,
                              arma::vec count, arma::vec M, arma::vec steps,
                              const unsigned int M_iter, const unsigned int cmin,
                              const unsigned int maxiter,
                              arma::vec start, const unsigned int burnin,
                              const bool progress = true){
  steps = steps - 1; // R --> C++ indexing
  vec inside;
  uvec inside_idx;
  mat mcmc;
  mat starts(steps.n_elem, A.n_cols);        // first row: A[1:steps[1],] satisfied
  for (unsigned int s = 0; s < steps.n_elem; s++){
    starts.row(s) = start_random(A, b, M(0), start).t(); // different starting values for each step
  }
  int i, from, iter = 0;  // from: can be negative!
  while(count.min() < cmin){
    Rcpp::checkUserInterrupt();
    if (progress && iter % (maxiter/100) == 0)
      Rcout << (iter== 0 ? " current cmin: " : " , ") << count.min();
    iter += 1;
    // find step with smallest count:
    i = count.index_min();

    // sampling for step with minimal counts
    if (i == 0){
      from = -1;   // i==0: start at A(0,:). note +1 below
      mcmc = rpdirichlet(M_iter, k + prior, options, true);
    }  else {
      from = steps(i-1);
      mcmc = sampling_mult(k, options, A.rows(0, from), b.subvec(0, from),
                            prior, M_iter, starts.row(i - 1).t(), burnin, false);
      starts.row(i - 1) = mcmc.row(M_iter - 1);
    }
    // count samples
    inside = inside_Ab(mcmc, A.rows(from + 1, steps(i)), b.subvec(from + 1, steps(i)));
    if (any(inside)){
      inside_idx = find(inside);
      starts.row(i) = mcmc.row(inside_idx(inside_idx.n_elem - 1));
    }
    count(i) += accu(inside);
    M(i) += M_iter;
    // compute precision
  }
  if (progress) Rcout << "\n";
  return results(count, M, steps + 1); // C++ --> R indexing
}



//#######################################################################
// BISECTION (in one dimension of a vector)

template <typename T>
double bisection(T f, NumericVector x, int i, double min, double max,
                 const double eps = 1e-10){

  // [robustness] check that lower/upper boundary are 0/1:
  x[i] = min;
  double f_min = as<double>(f(x)) - 0.50;
  x[i] = max;
  double f_max = as<double>(f(x)) - 0.50;
  if ( (f_min <= 0 && f_max <= 0) || (f_min >= 0 && f_max >= 0)){
    Rcout << "Bisection with respect to element [" << i+1 << "] on the interval [" << min << "," << max << "]\n";
    Rcout << "Current state of probability vector: " << x << "\n";
    stop("[Bisection algorithm]\n  Indicator function 'inside' does not have different values (0/1) for min/max.\n  Check whether inequality-constrained parameter space is convex!\n  (multiplicative constraints such as x[1]*x[2]<0.50 are in general not convex)");
  }

  while (min + eps < max) {
    double const mid = 0.5 * min + 0.5 * max;
    x[i] = mid;
    double const f_mid = as<double>(f(x)) - 0.50;

    if ((f_min < 0) == (f_mid < 0)) {
      min = mid;
      f_min = f_mid;
    } else {
      max = mid;
    }
  }

  return min;
}

// [[Rcpp::export]]
double bisection_r(Function f, NumericVector x, int i, double min, double max,
                   const double eps = 1e-10){
  return bisection<Function>(f, x, i, min, max, eps);
}

typedef SEXP (*funcPtr)(NumericVector);
// [[Rcpp::export]]
double bisection_cpp(SEXP f_, NumericVector x, int i, double min, double max,
                     const double eps = 1e-10){
  funcPtr f = *XPtr<funcPtr>(f_);
  return bisection<funcPtr>(f, x, i, min, max, eps);
}



//#######################################################################
// NONLINEAR CONSTRAINTS


typedef SEXP (*funcPtr)(NumericVector);
// [[Rcpp::export]]
NumericVector call_xptr(SEXP f_, NumericVector x){
  funcPtr f = *XPtr<funcPtr>(f_);
  NumericVector y = f(x);
  return y;
}

template <typename T>
arma::mat sampling_nonlin(const arma::vec& k, const arma::vec& options, T inside,
                          const arma::vec& prior, const unsigned int M, arma::vec start,
                          const unsigned int burnin = 5, const bool progress = true,
                          const double eps = 1e-10){
  unsigned int D = sum(options - 1);  // dimensions
  mat X(D, M + burnin);       // initialize posterior, column-major order
  X.col(0) = start;

  // last options for each multinomial (not provided by matrix A)
  uvec idx_J = conv_to<uvec>::from(cumsum(options) - 1);
  vec beta_J = rep_options(k(idx_J) + prior(idx_J), options - 1);
  vec beta_j = shed_options(k + prior, options);
  // lower/upper index within each multinomial:
  vec j_up  = rep_options(cumsum(options - 1) - 1, options - 1);
  vec j_low = j_up - rep_options(options - 2, options - 1);

  Progress p(M, progress);
  double bmax, bmin, s;
  uvec Apos, Aneg;
  IntegerVector idx = seq_len(D) - 1;
  unsigned int j = 0;
  for (unsigned int i = 1 ; i < M  + burnin; i++){
    p.increment();   // update progress bar
    if(i % 100 == 0) Rcpp::checkUserInterrupt();
    X.col(i) = X.col(i-1);
    idx = sample(idx, D, false);
    for (unsigned int m = 0; m < D; m++){
      j = idx(m);
      // scaling parameter of scaled truncated beta:
      s = 1 - accu(X.col(i).rows(j_low(j), j_up(j))) + X(j,i);
      // find truncation boundaries via bisection:
      bmin = 0.;
      bmax = s;

      X(j,i) = bmin;
      double check0 = as<double>(inside(wrap(X.col(i))));
      if (check0 != 1.){
        double rmin = bisection(inside, wrap(X.col(i)), j, bmin, X(j,i-1), eps);
        // approximation from below: add epsilon to ensure that inside(p)==TRUE
        bmin = fmax(0., rmin + eps);
      }

      X(j,i) = bmax;
      double check1 = as<double>(inside(wrap(X.col(i))));
      if (check1 != 1.){
        double rmax = bisection(inside, wrap(X.col(i)), j, X(j,i-1), bmax, eps);
        bmax = fmin(s, rmax);
      }
      X(j,i) = s * rbeta_trunc(beta_j(j), beta_J(j), bmin / s, bmax / s);
      double check_new = as<double>(inside(wrap(X.col(i))));
      if (check_new != 1.){
        Rcout << "index of updated parameter = " << j+1 ;
        Rcout << "\n  [scaling factor s = " << s << ";  lower bound = "<< bmin;
        Rcout <<";  upper bound = "<< bmax<<"]\n\ntheta[m]   = " << X.col(i).t();
        Rcout << "theta[m-1] = " << X.col(i-1).t() << "\n";
        stop("Gibbs sampler outside of truncated parameter space.\nPlease check whether restricted parameter space is convex.");
      }
    }
  }
  X.shed_cols(0,burnin - 1);
  return X.t();
}

// [[Rcpp::export]]
arma::mat sampling_nonlin_r(const arma::vec& k, const arma::vec& options, Function inside,
                            const arma::vec& prior, const unsigned int M, arma::vec start,
                            const unsigned int burnin = 5, const bool progress = true,
                            const double eps = 1e-10){
  return sampling_nonlin<Function>(k, options, inside, prior, M, start, burnin, progress, eps);
}

typedef SEXP (*funcPtr)(NumericVector);
// [[Rcpp::export]]
arma::mat sampling_nonlin_cpp(const arma::vec& k, const arma::vec& options, SEXP inside_,
                              const arma::vec& prior, const unsigned int M, arma::vec start,
                              const unsigned int burnin = 5, const bool progress = true,
                              const double eps = 1e-10){
  funcPtr inside = *XPtr<funcPtr>(inside_);
  return sampling_nonlin<funcPtr>(k, options, inside, prior, M, start, burnin, progress, eps);
}


typedef SEXP (*funcPtr)(NumericVector);
// [[Rcpp::export]]
NumericMatrix count_nonlin_cpp(const arma::vec& k, const arma::vec& options, SEXP inside_,
                               const arma::vec& prior, const unsigned int M,
                               const unsigned int batch, const bool progress = true){
  funcPtr inside = *XPtr<funcPtr>(inside_);

  Progress p(M/batch, progress);
  int count = 0, todo = M;
  mat X(batch, k.n_elem);
  while (todo > 0){
    p.increment();   // update progress bar
    Rcpp::checkUserInterrupt();
    // count prior and posterior samples that match constraints:
    unsigned int R = fmin(todo,batch);
    X = rpdirichlet(R, k + prior, options, true);
    for (unsigned int i = 0; i < R; i++){
      count = count  + as<double>(inside(wrap(X.row(i))));
    }
    todo = todo - batch;
  }
  return results(count, M, 1);
}
