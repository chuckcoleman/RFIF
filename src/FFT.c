/*
  FFT backend for RFIF.

  Preferred backend (when available): FFTW3 (fast).
  Fallback backend: self-contained FFT (portable, slower).

  The FFTW path is enabled when the package is configured with -DRFIF_USE_FFTW
  and linked with fftw3.
*/

#include "rfif_r.h"
#include <stdlib.h>
#include <math.h>
#include "FFT.h"

#ifdef RFIF_USE_FFTW
#include <fftw3.h>

double* realFFT(double *f, int N) {
  double *in = (double*)fftw_malloc(sizeof(double) * (size_t)N);
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (size_t)(N/2 + 1));
  double *re = (double*)malloc(sizeof(double) * (size_t)N);
  if (!in || !out || !re) {
    if (in) fftw_free(in);
    if (out) fftw_free(out);
    free(re);
    error("RFIF FFTW: allocation failed");
  }
  for (int i=0;i<N;i++) in[i] = f[i];

  fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
  if (!p) {
    fftw_free(in); fftw_free(out); free(re);
    error("RFIF FFTW: plan creation failed");
  }
  fftw_execute(p);

  re[0] = out[0][0];
  for (int k=1;k<N/2;k++) {
    re[k] = out[k][0];
    re[N-k] = out[k][0];
  }
  if (N % 2 == 0) re[N/2] = out[N/2][0];

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  return re;
}

fif_complex* fft_dir(double *f, int N) {
  double *in = (double*)fftw_malloc(sizeof(double) * (size_t)N);
  fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (size_t)(N/2 + 1));
  fif_complex *X = (fif_complex*)malloc(sizeof(fif_complex) * (size_t)N);
  if (!in || !out || !X) {
    if (in) fftw_free(in);
    if (out) fftw_free(out);
    free(X);
    error("RFIF FFTW: allocation failed");
  }
  for (int i=0;i<N;i++) in[i] = f[i];

  fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
  if (!p) {
    fftw_free(in); fftw_free(out); free(X);
    error("RFIF FFTW: plan creation failed");
  }
  fftw_execute(p);

  X[0].re = out[0][0]; X[0].im = out[0][1];
  for (int k=1;k<N/2;k++) {
    X[k].re = out[k][0]; X[k].im = out[k][1];
    X[N-k].re = out[k][0]; X[N-k].im = -out[k][1];
  }
  if (N % 2 == 0) {
    X[N/2].re = out[N/2][0]; X[N/2].im = out[N/2][1];
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  return X;
}

double* fft_inv(fif_complex *X, int N) {
  fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (size_t)(N/2 + 1));
  double *out = (double*)fftw_malloc(sizeof(double) * (size_t)N);
  double *f = (double*)malloc(sizeof(double) * (size_t)N);
  if (!in || !out || !f) {
    if (in) fftw_free(in);
    if (out) fftw_free(out);
    free(f);
    error("RFIF FFTW: allocation failed");
  }

  in[0][0] = X[0].re; in[0][1] = X[0].im;
  for (int k=1;k<N/2;k++) {
    in[k][0] = X[k].re; in[k][1] = X[k].im;
  }
  if (N % 2 == 0) {
    in[N/2][0] = X[N/2].re; in[N/2][1] = X[N/2].im;
  }

  fftw_plan p = fftw_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);
  if (!p) {
    fftw_free(in); fftw_free(out); free(f);
    error("RFIF FFTW: plan creation failed");
  }
  fftw_execute(p);

  for (int i=0;i<N;i++) f[i] = out[i] / (double)N;

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  return f;
}

#else  /* RFIF_USE_FFTW */

/* ---- Fallback: self-contained FFT (portable, slower) ---- */
typedef struct { double re; double im; } cpx;

static inline cpx c_add(cpx a, cpx b){ cpx r = {a.re+b.re, a.im+b.im}; return r; }
static inline cpx c_sub(cpx a, cpx b){ cpx r = {a.re-b.re, a.im-b.im}; return r; }
static inline cpx c_mul(cpx a, cpx b){
  cpx r = {a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re};
  return r;
}
static inline cpx c_exp_i(double theta){
  cpx r = {cos(theta), sin(theta)};
  return r;
}

/* Return smallest factor >1 of n, or n if prime. */
static int smallest_factor(int n){
  if(n % 2 == 0) return 2;
  for(int p=3; p*(long)p <= n; p += 2){
    if(n % p == 0) return p;
  }
  return n;
}

/* Naive DFT (used for prime sizes). inverse=0 forward, inverse=1 inverse (unnormalized). */
static void dft_naive(const cpx *in, cpx *out, int n, int inverse){
  double sign = inverse ? 1.0 : -1.0;
  for(int k=0;k<n;k++){
    double sum_re=0.0, sum_im=0.0;
    for(int t=0;t<n;t++){
      double ang = sign * 2.0 * M_PI * (double)k * (double)t / (double)n;
      double ca = cos(ang), sa = sin(ang);
      /* in[t] * exp(i ang) */
      sum_re += in[t].re*ca - in[t].im*sa;
      sum_im += in[t].re*sa + in[t].im*ca;
    }
    out[k].re = sum_re;
    out[k].im = sum_im;
  }
}

/* Recursive mixed-radix FFT. */
static void fft_rec(const cpx *in, cpx *out, int n, int inverse){
  if(n == 1){
    out[0] = in[0];
    return;
  }
  int p = smallest_factor(n);
  if(p == n){
    dft_naive(in, out, n, inverse);
    return;
  }
  int m = n / p;

  /* Step 1: compute p FFTs of size m on decimated sequences */
  cpx *tmp = (cpx*)malloc(sizeof(cpx) * (size_t)n);
  cpx *buf_in  = (cpx*)malloc(sizeof(cpx) * (size_t)m);
  cpx *buf_out = (cpx*)malloc(sizeof(cpx) * (size_t)m);
  if(!tmp || !buf_in || !buf_out){
    free(tmp); free(buf_in); free(buf_out);
    error("RFIF FFT: out of memory");
  }

  for(int r=0;r<p;r++){
    for(int j=0;j<m;j++){
      buf_in[j] = in[j*p + r];
    }
    fft_rec(buf_in, buf_out, m, inverse);
    for(int j=0;j<m;j++){
      tmp[r*m + j] = buf_out[j];
    }
  }

  /* Step 2: combine with twiddle factors */
  double sign = inverse ? 1.0 : -1.0;
  for(int k=0;k<n;k++){
    int q = k / m;      /* 0..p-1 */
    int s = k % m;      /* 0..m-1 */
    cpx acc = {0.0, 0.0};
    for(int r=0;r<p;r++){
      /* twiddle: exp(sign * 2pi i * r * k / n) */
      double ang = sign * 2.0 * M_PI * (double)r * (double)k / (double)n;
      cpx w = c_exp_i(ang);
      acc = c_add(acc, c_mul(tmp[r*m + s], w));
    }
    out[k] = acc;
  }

  free(tmp);
  free(buf_in);
  free(buf_out);
}

static void rfif_fft_inplace(double *re, double *im, int n, int inverse){
  cpx *in  = (cpx*)malloc(sizeof(cpx) * (size_t)n);
  cpx *out = (cpx*)malloc(sizeof(cpx) * (size_t)n);
  if(!in || !out){
    free(in); free(out);
    error("RFIF FFT: allocation failed");
  }
  for(int i=0;i<n;i++){
    in[i].re = re[i];
    in[i].im = im ? im[i] : 0.0;
  }

  fft_rec(in, out, n, inverse);

  /* Normalize inverse */
  if(inverse){
    for(int i=0;i<n;i++){
      re[i] = out[i].re / (double)n;
      if(im) im[i] = out[i].im / (double)n;
    }
  }else{
    for(int i=0;i<n;i++){
      re[i] = out[i].re;
      if(im) im[i] = out[i].im;
    }
  }

  free(in);
  free(out);
}

double* realFFT(double *f, int N) {
  double *re = (double*) malloc(sizeof(double) * (size_t)N);
  double *im = (double*) calloc((size_t)N, sizeof(double));
  if (!re || !im) {
    free(re); free(im);
    error("Allocation failed in realFFT");
  }
  for (int i = 0; i < N; ++i) re[i] = f[i];

  rfif_fft_inplace(re, im, N, 0);

  free(im);
  return re;
}

fif_complex* fft_dir(double *f, int N) {
  double *re = (double*) malloc(sizeof(double) * (size_t)N);
  double *im = (double*) calloc((size_t)N, sizeof(double));
  fif_complex *out = (fif_complex*) malloc(sizeof(fif_complex) * (size_t)N);

  if (!re || !im || !out) {
    free(re); free(im); free(out);
    error("Allocation failed in fft_dir");
  }

  for (int i = 0; i < N; ++i) re[i] = f[i];

  rfif_fft_inplace(re, im, N, 0);

  for (int i = 0; i < N; ++i) {
    out[i].re = re[i];
    out[i].im = im[i];
  }

  free(re);
  free(im);
  return out;
}

double* fft_inv(fif_complex *X, int N) {
  double *re = (double*) malloc(sizeof(double) * (size_t)N);
  double *im = (double*) malloc(sizeof(double) * (size_t)N);
  double *out = (double*) malloc(sizeof(double) * (size_t)N);

  if (!re || !im || !out) {
    free(re); free(im); free(out);
    error("Allocation failed in fft_inv");
  }

  for (int i = 0; i < N; ++i) {
    re[i] = X[i].re;
    im[i] = X[i].im;
  }

  rfif_fft_inplace(re, im, N, 1);

  for (int i = 0; i < N; ++i) out[i] = re[i];

  free(re);
  free(im);
  return out;
}

#endif /* RFIF_USE_FFTW */
