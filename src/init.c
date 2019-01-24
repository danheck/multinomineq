#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _multinomineq_adj_iterative(SEXP, SEXP, SEXP);
extern SEXP _multinomineq_bisection_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_bisection_r(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_call_xptr(SEXP, SEXP);
extern SEXP _multinomineq_count_auto_bin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_count_auto_mult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_count_bin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_count_mult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_count_nonlin_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_count_samples(SEXP, SEXP, SEXP);
extern SEXP _multinomineq_count_stepwise_bin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_count_stepwise_multi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_inside_Ab(SEXP, SEXP, SEXP);
extern SEXP _multinomineq_ppp_bin(SEXP, SEXP, SEXP);
extern SEXP _multinomineq_ppp_mult(SEXP, SEXP, SEXP);
extern SEXP _multinomineq_rbeta_trunc(SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_rdirichlet(SEXP, SEXP);
extern SEXP _multinomineq_rep_options(SEXP, SEXP);
extern SEXP _multinomineq_rgamma_trunc(SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_rpdirichlet(SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_rpm_mat(SEXP, SEXP, SEXP);
extern SEXP _multinomineq_sampling_bin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_sampling_hitandrun(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_sampling_mult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_sampling_nonlin_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_sampling_nonlin_r(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_start_random(SEXP, SEXP, SEXP, SEXP);
extern SEXP _multinomineq_sum_options(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_multinomineq_adj_iterative",        (DL_FUNC) &_multinomineq_adj_iterative,         3},
    {"_multinomineq_bisection_cpp",        (DL_FUNC) &_multinomineq_bisection_cpp,         6},
    {"_multinomineq_bisection_r",          (DL_FUNC) &_multinomineq_bisection_r,           6},
    {"_multinomineq_call_xptr",            (DL_FUNC) &_multinomineq_call_xptr,             2},
    {"_multinomineq_count_auto_bin",       (DL_FUNC) &_multinomineq_count_auto_bin,       14},
    {"_multinomineq_count_auto_mult",      (DL_FUNC) &_multinomineq_count_auto_mult,      14},
    {"_multinomineq_count_bin",            (DL_FUNC) &_multinomineq_count_bin,             8},
    {"_multinomineq_count_mult",           (DL_FUNC) &_multinomineq_count_mult,            8},
    {"_multinomineq_count_nonlin_cpp",     (DL_FUNC) &_multinomineq_count_nonlin_cpp,      7},
    {"_multinomineq_count_samples",        (DL_FUNC) &_multinomineq_count_samples,         3},
    {"_multinomineq_count_stepwise_bin",   (DL_FUNC) &_multinomineq_count_stepwise_bin,   11},
    {"_multinomineq_count_stepwise_multi", (DL_FUNC) &_multinomineq_count_stepwise_multi, 11},
    {"_multinomineq_inside_Ab",            (DL_FUNC) &_multinomineq_inside_Ab,             3},
    {"_multinomineq_ppp_bin",              (DL_FUNC) &_multinomineq_ppp_bin,               3},
    {"_multinomineq_ppp_mult",             (DL_FUNC) &_multinomineq_ppp_mult,              3},
    {"_multinomineq_rbeta_trunc",          (DL_FUNC) &_multinomineq_rbeta_trunc,           4},
    {"_multinomineq_rdirichlet",           (DL_FUNC) &_multinomineq_rdirichlet,            2},
    {"_multinomineq_rep_options",          (DL_FUNC) &_multinomineq_rep_options,           2},
    {"_multinomineq_rgamma_trunc",         (DL_FUNC) &_multinomineq_rgamma_trunc,          4},
    {"_multinomineq_rpdirichlet",          (DL_FUNC) &_multinomineq_rpdirichlet,           4},
    {"_multinomineq_rpm_mat",              (DL_FUNC) &_multinomineq_rpm_mat,               3},
    {"_multinomineq_sampling_bin",         (DL_FUNC) &_multinomineq_sampling_bin,          9},
    {"_multinomineq_sampling_hitandrun",   (DL_FUNC) &_multinomineq_sampling_hitandrun,    6},
    {"_multinomineq_sampling_mult",        (DL_FUNC) &_multinomineq_sampling_mult,         9},
    {"_multinomineq_sampling_nonlin_cpp",  (DL_FUNC) &_multinomineq_sampling_nonlin_cpp,   9},
    {"_multinomineq_sampling_nonlin_r",    (DL_FUNC) &_multinomineq_sampling_nonlin_r,     9},
    {"_multinomineq_start_random",         (DL_FUNC) &_multinomineq_start_random,          4},
    {"_multinomineq_sum_options",          (DL_FUNC) &_multinomineq_sum_options,           2},
    {NULL, NULL, 0}
};

void R_init_multinomineq(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
