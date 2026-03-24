#ifndef RFIF_FFT_H
#define RFIF_FFT_H

typedef struct fif_complex {
    double re;
    double im;
} fif_complex;

double* realFFT(double *f, int N);
fif_complex* fft_dir(double *f, int N);
double* fft_inv(fif_complex *X, int N);

#endif
