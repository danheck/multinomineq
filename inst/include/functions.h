using namespace Rcpp;
using namespace arma;

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

// functions that are used in multiple scripts:

NumericMatrix results(const vec& count, const vec& M, const vec& steps);
NumericMatrix results(const unsigned int count, const unsigned int M, const unsigned int steps);

double x2(const vec& o, const vec& e);

vec inside_Ab(const mat& X, const mat& A, const vec& b);
int count_samples(const mat& X, const mat& A, const vec& b);
int count_samples(const vec& X, const mat& A, const vec& b);

vec sort_steps(vec steps, const unsigned int max);

double rbeta_trunc(const double shape1, const double shape2,
                   const double min, const double max);

vec start_random(const mat& A, const vec& b, const unsigned int M, vec start);

#endif

