#include "rfif_r.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "Fif.h"

// Free FIF output list produced by FIF_v2_1.
// The head struct is returned by value, so we free head.dati but not the head itself.
// Subsequent nodes are heap-allocated.
static void free_fif_output(Fif_t *head) {
  if (!head) return;

  // free head data
  if (head->dati) {
    free(head->dati);
    head->dati = NULL;
  }

  Fif_t *cur = head->next;
  head->next = NULL;

  while (cur) {
    Fif_t *nxt = cur->next;
    if (cur->dati) free(cur->dati);
    free(cur);
    cur = nxt;
  }
}

SEXP r_fifc_run(SEXP x) {
  if (!Rf_isReal(x)) Rf_error("x must be a numeric vector");
  R_xlen_t N_x = XLENGTH(x);
  if (N_x < 2) Rf_error("x too short");
  if (N_x > INT_MAX) Rf_error("x too long for this implementation");

  int N = (int)N_x;
  const double *f_in = REAL(x);

  // CRAN-hardening: pass a private copy into the C core in case it mutates input
  double *f_work = (double*)malloc((size_t)N * sizeof(double));
  if (!f_work) Rf_error("RFIF: allocation failed (f_work)");
  memcpy(f_work, f_in, (size_t)N * sizeof(double));

  int nimf = 0;
  Fif_t out = FIF_v2_1(f_work, N, &nimf);

  // f_work is no longer needed beyond this point
  free(f_work);
  f_work = NULL;

  if (nimf <= 0 || !out.dati) {
    // No oscillatory modes found: return empty IMF matrix and treat input as residual.
    free_fif_output(&out);

    SEXP imfs = PROTECT(Rf_allocMatrix(REALSXP, 0, N));
    SEXP residual = PROTECT(Rf_allocVector(REALSXP, N));
    double *rres = REAL(residual);
    const double *xin = REAL(x);
    for (int i = 0; i < N; ++i) rres[i] = xin[i];

    SEXP nimfS = PROTECT(Rf_ScalarInteger(0));

    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, imfs);
    SET_VECTOR_ELT(ans, 1, residual);
    SET_VECTOR_ELT(ans, 2, nimfS);

    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, Rf_mkChar("imfs"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("residual"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("nimf"));
    Rf_setAttrib(ans, R_NamesSymbol, nms);

    UNPROTECT(6);
    return ans;
  }

  // Allocate IMFs matrix [nimf x N]
  SEXP imfs = PROTECT(Rf_allocMatrix(REALSXP, nimf, N));
  double *rimfs = REAL(imfs);

  Fif_t *cur = &out;
  for (int k = 0; k < nimf; ++k) {
    if (!cur || !cur->dati) {
      free_fif_output(&out);
      UNPROTECT(1);
      Rf_error("FIF output list shorter than nimf (k=%d)", k);
    }
    for (int i = 0; i < N; ++i) {
      rimfs[k + nimf * i] = cur->dati[i];  // column-major
    }
    cur = cur->next;
  }

  // residual = x - sum(IMFs)
  SEXP resid = PROTECT(Rf_allocVector(REALSXP, N));
  double *rres = REAL(resid);
  for (int i = 0; i < N; ++i) {
    double s = 0.0;
    for (int k = 0; k < nimf; ++k) s += rimfs[k + nimf * i];
    rres[i] = f_in[i] - s;
  }

  SEXP outlist = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(outlist, 0, imfs);
  SET_VECTOR_ELT(outlist, 1, resid);
  SET_VECTOR_ELT(outlist, 2, Rf_ScalarInteger(nimf));

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("imfs"));
  SET_STRING_ELT(names, 1, Rf_mkChar("residual"));
  SET_STRING_ELT(names, 2, Rf_mkChar("nimf"));
  Rf_setAttrib(outlist, R_NamesSymbol, names);

  // free native output list (we've copied into R)
  free_fif_output(&out);

  UNPROTECT(4);
  return outlist;
}
