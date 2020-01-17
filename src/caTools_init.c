#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void cumsum_exact(void *, void *, void *);
extern void imwritegif(void *, void *, void *, void *, void *);
extern void runmad(void *, void *, void *, void *, void *);
extern void runmax(void *, void *, void *, void *);
extern void runmean(void *, void *, void *, void *);
extern void runmean_exact(void *, void *, void *, void *);
extern void runmean_lite(void *, void *, void *, void *);
extern void runmin(void *, void *, void *, void *);
extern void runquantile(void *, void *, void *, void *, void *, void *, void *);
extern void runsd(void *, void *, void *, void *, void *);
extern void sum_exact(void *, void *, void *);

/* .Call calls */
extern SEXP imreadgif(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"cumsum_exact",  (DL_FUNC) &cumsum_exact,  3},
    {"imwritegif",    (DL_FUNC) &imwritegif,    5},
    {"runmad",        (DL_FUNC) &runmad,        5},
    {"runmax",        (DL_FUNC) &runmax,        4},
    {"runmean",       (DL_FUNC) &runmean,       4},
    {"runmean_exact", (DL_FUNC) &runmean_exact, 4},
    {"runmean_lite",  (DL_FUNC) &runmean_lite,  4},
    {"runmin",        (DL_FUNC) &runmin,        4},
    {"runquantile",   (DL_FUNC) &runquantile,   7},
    {"runsd",         (DL_FUNC) &runsd,         5},
    {"sum_exact",     (DL_FUNC) &sum_exact,     3},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"imreadgif", (DL_FUNC) &imreadgif, 3},
    {NULL, NULL, 0}
};

void R_init_caTools(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
