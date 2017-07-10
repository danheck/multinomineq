#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP stratsel_adj_iterative(SEXP, SEXP, SEXP);
extern SEXP stratsel_count_samples(SEXP, SEXP, SEXP);
extern SEXP stratsel_encompassing_bf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_encompassing_stepwise(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stratsel_sampling_binary_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"stratsel_adj_iterative",         (DL_FUNC) &stratsel_adj_iterative,         3},
    {"stratsel_count_samples",         (DL_FUNC) &stratsel_count_samples,         3},
    {"stratsel_encompassing_bf",       (DL_FUNC) &stratsel_encompassing_bf,       7},
    {"stratsel_encompassing_stepwise", (DL_FUNC) &stratsel_encompassing_stepwise, 8},
    {"stratsel_sampling_binary_cpp",   (DL_FUNC) &stratsel_sampling_binary_cpp,   7},
    {NULL, NULL, 0}
};

void R_init_stratsel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
