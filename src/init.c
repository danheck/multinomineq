#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP stratsel_adj_iterative(SEXP, SEXP, SEXP);
extern SEXP stratsel_count_binomial_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_count_multinomial_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_count_samples(SEXP, SEXP, SEXP);
extern SEXP stratsel_count_stepwise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_rbeta_trunc(SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_rdirichlet(SEXP, SEXP);
extern SEXP stratsel_rep_options(SEXP, SEXP);
extern SEXP stratsel_rpdirichlet(SEXP, SEXP, SEXP);
extern SEXP stratsel_rpdirichlet_free(SEXP, SEXP, SEXP);
extern SEXP stratsel_sampling_binomial_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_sampling_multinomial_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_shed_options(SEXP, SEXP);
extern SEXP stratsel_start_random(SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_sum_options(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"stratsel_adj_iterative",            (DL_FUNC) &stratsel_adj_iterative,            3},
    {"stratsel_count_binomial_cpp",       (DL_FUNC) &stratsel_count_binomial_cpp,       7},
    {"stratsel_count_multinomial_cpp",    (DL_FUNC) &stratsel_count_multinomial_cpp,    7},
    {"stratsel_count_samples",            (DL_FUNC) &stratsel_count_samples,            3},
    {"stratsel_count_stepwise",           (DL_FUNC) &stratsel_count_stepwise,           8},
    {"stratsel_rbeta_trunc",              (DL_FUNC) &stratsel_rbeta_trunc,              4},
    {"stratsel_rdirichlet",               (DL_FUNC) &stratsel_rdirichlet,               2},
    {"stratsel_rep_options",              (DL_FUNC) &stratsel_rep_options,              2},
    {"stratsel_rpdirichlet",              (DL_FUNC) &stratsel_rpdirichlet,              3},
    {"stratsel_rpdirichlet_free",         (DL_FUNC) &stratsel_rpdirichlet_free,         3},
    {"stratsel_sampling_binomial_cpp",    (DL_FUNC) &stratsel_sampling_binomial_cpp,    8},
    {"stratsel_sampling_multinomial_cpp", (DL_FUNC) &stratsel_sampling_multinomial_cpp, 8},
    {"stratsel_shed_options",             (DL_FUNC) &stratsel_shed_options,             2},
    {"stratsel_start_random",             (DL_FUNC) &stratsel_start_random,             4},
    {"stratsel_sum_options",              (DL_FUNC) &stratsel_sum_options,              2},
    {NULL, NULL, 0}
};

void R_init_stratsel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
