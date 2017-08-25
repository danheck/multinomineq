#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _stratsel_adj_iterative(SEXP, SEXP, SEXP);
extern SEXP _stratsel_count_binomial_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_count_multinomial_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_count_samples(SEXP, SEXP, SEXP);
extern SEXP _stratsel_count_step(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_count_stepwise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_count_stepwise_multi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_inside_Ab(SEXP, SEXP, SEXP);
extern SEXP _stratsel_ppp_bin(SEXP, SEXP, SEXP);
extern SEXP _stratsel_ppp_mult(SEXP, SEXP, SEXP);
extern SEXP _stratsel_rbeta_trunc(SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_rdirichlet(SEXP, SEXP);
extern SEXP _stratsel_rep_options(SEXP, SEXP);
extern SEXP _stratsel_rpdirichlet(SEXP, SEXP, SEXP);
extern SEXP _stratsel_rpdirichlet_free(SEXP, SEXP, SEXP);
extern SEXP _stratsel_rpm_mat(SEXP, SEXP, SEXP);
extern SEXP _stratsel_sampling_binomial_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_sampling_hitandrun(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_sampling_multinomial_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_shed_options(SEXP, SEXP);
extern SEXP _stratsel_start_random(SEXP, SEXP, SEXP, SEXP);
extern SEXP _stratsel_sum_options(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_stratsel_adj_iterative",            (DL_FUNC) &_stratsel_adj_iterative,             3},
    {"_stratsel_count_binomial_cpp",       (DL_FUNC) &_stratsel_count_binomial_cpp,        8},
    {"_stratsel_count_multinomial_cpp",    (DL_FUNC) &_stratsel_count_multinomial_cpp,     8},
    {"_stratsel_count_samples",            (DL_FUNC) &_stratsel_count_samples,             3},
    {"_stratsel_count_step",               (DL_FUNC) &_stratsel_count_step,               10},
    {"_stratsel_count_stepwise",           (DL_FUNC) &_stratsel_count_stepwise,           10},
    {"_stratsel_count_stepwise_multi",     (DL_FUNC) &_stratsel_count_stepwise_multi,     10},
    {"_stratsel_inside_Ab",                (DL_FUNC) &_stratsel_inside_Ab,                 3},
    {"_stratsel_ppp_bin",                  (DL_FUNC) &_stratsel_ppp_bin,                   3},
    {"_stratsel_ppp_mult",                 (DL_FUNC) &_stratsel_ppp_mult,                  3},
    {"_stratsel_rbeta_trunc",              (DL_FUNC) &_stratsel_rbeta_trunc,               4},
    {"_stratsel_rdirichlet",               (DL_FUNC) &_stratsel_rdirichlet,                2},
    {"_stratsel_rep_options",              (DL_FUNC) &_stratsel_rep_options,               2},
    {"_stratsel_rpdirichlet",              (DL_FUNC) &_stratsel_rpdirichlet,               3},
    {"_stratsel_rpdirichlet_free",         (DL_FUNC) &_stratsel_rpdirichlet_free,          3},
    {"_stratsel_rpm_mat",                  (DL_FUNC) &_stratsel_rpm_mat,                   3},
    {"_stratsel_sampling_binomial_cpp",    (DL_FUNC) &_stratsel_sampling_binomial_cpp,     9},
    {"_stratsel_sampling_hitandrun",       (DL_FUNC) &_stratsel_sampling_hitandrun,        6},
    {"_stratsel_sampling_multinomial_cpp", (DL_FUNC) &_stratsel_sampling_multinomial_cpp,  9},
    {"_stratsel_shed_options",             (DL_FUNC) &_stratsel_shed_options,              2},
    {"_stratsel_start_random",             (DL_FUNC) &_stratsel_start_random,              4},
    {"_stratsel_sum_options",              (DL_FUNC) &_stratsel_sum_options,               2},
    {NULL, NULL, 0}
};

void R_init_stratsel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
