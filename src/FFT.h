#ifndef RFIF_FFT_H
#define RFIF_FFT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct fif_complex {
    double re;
    double im;
} fif_complex;

double *realFFT(double *x, int N);
fif_complex *fft_dir(double *x, int N);
double *fft_inv(fif_complex *X, int N);

#ifdef __cplusplus
}
#endif

#endif