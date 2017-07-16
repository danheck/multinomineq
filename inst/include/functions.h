using namespace Rcpp;
using namespace arma;

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

// functions that are used in multiple scripts:

int count_samples(mat X, mat A, vec b);
int count_samples(vec X, mat A, vec b);

double rbeta_trunc(double shape1, double shape2, double min, double max);

vec start_random(mat A, vec b, int M, vec start);

#endif
