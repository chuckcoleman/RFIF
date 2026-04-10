#ifndef RFIF_FFT_H
#define RFIF_FFT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct fif_complex {
    double re;
    double im;
} fif_complex;

/* Real part of forward transform of a real-valued signal.
   Returns an allocated vector of length N, or NULL on failure. */
double *realFFT(double *x, int N);

/* Forward transform of a real-valued signal into complex spectrum.
   Returns an allocated vector of length N, or NULL on failure. */
fif_complex *fft_dir(double *x, int N);

/* Inverse transform from complex spectrum to real-valued signal.
   Returns an allocated vector of length N, or NULL on failure. */
double *fft_inv(fif_complex *X, int N);

#ifdef __cplusplus
}
#endif

#endif