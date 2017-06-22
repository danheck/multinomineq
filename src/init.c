#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP stratsel_adj_iterative(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"stratsel_adj_iterative", (DL_FUNC) &stratsel_adj_iterative, 3},
    {NULL, NULL, 0}
};

void R_init_stratsel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
