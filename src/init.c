#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void pgmm_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"pgmm_c", (DL_FUNC) &pgmm_c, 13},
    {NULL, NULL, 0}
};

void R_init_pgmm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
