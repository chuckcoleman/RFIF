/* 
    FFT backend for RFIF

    Provides a single wrapper API for:
      - FFTW-backed execution when RFIF_USE_FFTW is defined
      - portable direct-DFT fallback otherwise

    Public API expected by FFT.h:
      double*      realFFT(double *x, int N);
      fif_complex* fft_dir(double *x, int N);
      double*      fft_inv(fif_complex *X, int N);
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FFT.h"

#ifdef RFIF_USE_FFTW
#include <fftw3.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ----------------------------- */
/* Internal helpers (fallback)   */
/* ----------------------------- */

/* ----------------------------- */
/* realFFT                       */
/* ----------------------------- */
/*
   Returns the real part of the forward transform of x.

   Caller owns the returned buffer and must free() it.
*/
double *realFFT(double *x, int N)
{
    double *out = NULL;

    if (x == NULL || N <= 0)
        return NULL;

    out = (double *)calloc((size_t)N, sizeof(double));
    if (out == NULL)
        return NULL;

#ifdef RFIF_USE_FFTW
    {
        fftw_complex *in = NULL, *freq = NULL;
        fftw_plan plan;

        in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (size_t)N);
        freq = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (size_t)N);
        if (in == NULL || freq == NULL)
        {
            if (in != NULL) fftw_free(in);
            if (freq != NULL) fftw_free(freq);
            free(out);
            return NULL;
        }

        for (int k = 0; k < N; ++k)
        {
            in[k][0] = x[k];
            in[k][1] = 0.0;
        }

        plan = fftw_plan_dft_1d(N, in, freq, FFTW_FORWARD, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            fftw_free(in);
            fftw_free(freq);
            free(out);
            return NULL;
        }

        fftw_execute(plan);

        for (int k = 0; k < N; ++k)
            out[k] = freq[k][0];

        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(freq);
    }
#else
    {
        for (int k = 0; k < N; ++k)
        {
            double sum_re = 0.0;

            for (int n = 0; n < N; ++n)
            {
                double angle = -2.0 * M_PI * (double)k * (double)n / (double)N;
                sum_re += x[n] * cos(angle);
            }

            out[k] = sum_re;
        }
    }
#endif

    return out;
}

/* ----------------------------- */
/* fft_dir                       */
/* ----------------------------- */
/*
   Forward transform of a real-valued signal into fif_complex.

   Caller owns the returned buffer and must free() it.
*/
fif_complex *fft_dir(double *x, int N)
{
    fif_complex *out = NULL;

    if (x == NULL || N <= 0)
        return NULL;

    out = (fif_complex *)calloc((size_t)N, sizeof(fif_complex));
    if (out == NULL)
        return NULL;

#ifdef RFIF_USE_FFTW
    {
        fftw_complex *in = NULL, *freq = NULL;
        fftw_plan plan;

        in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (size_t)N);
        freq = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (size_t)N);
        if (in == NULL || freq == NULL)
        {
            if (in != NULL) fftw_free(in);
            if (freq != NULL) fftw_free(freq);
            free(out);
            return NULL;
        }

        for (int n = 0; n < N; ++n)
        {
            in[n][0] = x[n];
            in[n][1] = 0.0;
        }

        plan = fftw_plan_dft_1d(N, in, freq, FFTW_FORWARD, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            fftw_free(in);
            fftw_free(freq);
            free(out);
            return NULL;
        }

        fftw_execute(plan);

        for (int k = 0; k < N; ++k)
        {
            out[k].re = freq[k][0];
            out[k].im = freq[k][1];
        }

        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(freq);
    }
#else
    {
        for (int k = 0; k < N; ++k)
        {
            double sum_re = 0.0;
            double sum_im = 0.0;

            for (int n = 0; n < N; ++n)
            {
                double angle = -2.0 * M_PI * (double)k * (double)n / (double)N;
                double ca = cos(angle);
                double sa = sin(angle);

                sum_re += x[n] * ca;
                sum_im += x[n] * sa;
            }

            out[k].re = sum_re;
            out[k].im = sum_im;
        }
    }
#endif

    return out;
}

/* ----------------------------- */
/* fft_inv                       */
/* ----------------------------- */
/*
   Inverse transform from fif_complex to a real-valued signal.

   Caller owns the returned buffer and must free() it.
*/
double *fft_inv(fif_complex *X, int N)
{
    double *out = NULL;

    if (X == NULL || N <= 0)
        return NULL;

    out = (double *)calloc((size_t)N, sizeof(double));
    if (out == NULL)
        return NULL;

#ifdef RFIF_USE_FFTW
    {
        fftw_complex *in = NULL, *time = NULL;
        fftw_plan plan;

        in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (size_t)N);
        time = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (size_t)N);
        if (in == NULL || time == NULL)
        {
            if (in != NULL) fftw_free(in);
            if (time != NULL) fftw_free(time);
            free(out);
            return NULL;
        }

        for (int k = 0; k < N; ++k)
        {
            in[k][0] = X[k].re;
            in[k][1] = X[k].im;
        }

        plan = fftw_plan_dft_1d(N, in, time, FFTW_BACKWARD, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            fftw_free(in);
            fftw_free(time);
            free(out);
            return NULL;
        }

        fftw_execute(plan);

        for (int n = 0; n < N; ++n)
            out[n] = time[n][0] / (double)N;

        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(time);
    }
#else
    {
        for (int n = 0; n < N; ++n)
        {
            double sum_re = 0.0;

            for (int k = 0; k < N; ++k)
            {
                double angle = 2.0 * M_PI * (double)k * (double)n / (double)N;
                double ca = cos(angle);
                double sa = sin(angle);

                sum_re += X[k].re * ca - X[k].im * sa;
            }

            out[n] = sum_re / (double)N;
        }
    }
#endif

    return out;
}