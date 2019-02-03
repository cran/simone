#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ARLasso(void *, void *, void *, void *, void *, void *);
extern void GLasso(void *, void *, void *, void *, void *, void *, void *, void *);
extern void Lasso(void *, void *, void *, void *, void *);
extern void vertex_coord_C(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"ARLasso",        (DL_FUNC) &ARLasso,        6},
    {"GLasso",         (DL_FUNC) &GLasso,         8},
    {"Lasso",          (DL_FUNC) &Lasso,          5},
    {"vertex_coord_C", (DL_FUNC) &vertex_coord_C, 9},
    {NULL, NULL, 0}
};

void R_init_simone(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
