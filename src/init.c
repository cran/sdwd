#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(sdwdnet)(double *lam2, int *nobs, int *nvars, double *x, double *y, int *jd, 
    double *pf, double *pf2, int *dfmax, int *pmax, int *nlam, double *flmin, double *ulam, 
    double *eps, int *isd, int *maxit, int *strong, int *nalam, double *b0, double *beta, 
    int *ibeta, int *nbeta, double *alam, int *npass, int *jerr);
extern void F77_NAME(sdwdnetpath)(double *lam2, double *maj, int *nobs, int *nvars, double *x, 
    double *y, int *ju, double *pf, double *pf2, int *dfmax, int *pmax, int *nlam, double *flmin, 
    double *ulam, double *eps, int *maxit, int *strong, int *nalam, double *b0, double *beta, 
    int *m, int *nbeta, double *alam, int *npass, int *jerr);
extern void F77_NAME(standard)(int *nobs, int *nvars, double *x, int *ju, int *isd, double *xmean, 
    double *xnorm, double *maj);
extern void F77_NAME(chkvars)(int *nobs, int *nvars, double *x, int *ju);
extern void F77_NAME(dwddrv)(int *nobs, int *nvars, double *x, double *y, double *r, double *vl);

static const R_FortranMethodDef FortranEntries[] = {
    {"sdwdnet",      (DL_FUNC) &F77_SUB(sdwdnet),  25},
    {"sdwdnetpath",  (DL_FUNC) &F77_SUB(sdwdnetpath),  25},
    {"standard",     (DL_FUNC) &F77_SUB(standard), 8},
    {"chkvars",      (DL_FUNC) &F77_SUB(chkvars),  4},
    {"dwddrv",       (DL_FUNC) &F77_SUB(dwddrv),  6},
    {NULL, NULL, 0}
};

void R_init_sdwd(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
