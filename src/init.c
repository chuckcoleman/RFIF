#include "rfif_r.h"
#include <R_ext/Rdynload.h>

SEXP r_fifc_run(SEXP x);

static const R_CallMethodDef CallEntries[] = {
  {"r_fifc_run", (DL_FUNC) &r_fifc_run, 1},
  {NULL, NULL, 0}
};

void R_init_RFIF(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
