/* 
    C implementation of Fast Iterative Filtering (FIF)
    
    written based on  FIF_v2_13.m (MATLAB VERSION)

    Authors: Igor Bertello, Emanuele Papini, Antonio Cicone
    Affiliation(s): IAPS - INAF, University of L'Aquila (Italy)
    
    Revised:  Chuck Coleman
    Affiliation:  Timely Analytics, LLC

    Dependencies: FFTW3
*/

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <fftw3.h>
#include "FFT.h"
#include "Fif.h"
#include "interp.h"
#include <stdbool.h>

#define DELTA 0.001
#define MaxInner 200
#define Xi 1.6 
#define ExtPoints 3

// the alpha value is set by default to average value (see line 490)
int comp(const void *elem1, const void *elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

Maxmins MaxMin;

double *fappo; // vettore di appoggio

//*********************** getMask *******************************
double *getMask(double *y, int n, double k, int *dim_a, double tol)
{
    /*
% Rescale the mask y so that its length becomes 2*k+1.
% k could be an integer or not an integer.
% y is the area under the curve for each bar
*/
    int lim = 0;
    double m = 0.0;
    double *a = NULL, *f_appo = NULL;
    double s = 0.0, s2 = 0.0, t = 0.0, t1 = 0.0, t2 = 0.0, c = 0.0;
    double parz = 0.0, extra = 0.0, new_k = 0.0, dx = 0.0, dy = 0.0;
    double Norm1 = 0.0;
    bool f_appoB = 0;

    if (y == NULL || dim_a == NULL || n <= 0) {
        return 0;
    }

    m = (double)(n - 1) / 2.;

    if (k <= m)
    {
        if ((int)k == k)
        {
            lim = 2 * (int)k + 1;

            a = (double *)calloc((size_t)lim, sizeof(double));
            if (a == NULL)
                return 0;

            for (int i = 1; i < lim + 1; i++)
            {
                int is, it, jstart, jend_excl;
                s = ((double)(i - 1) * (2 * m + 1)) / (double)(lim) + 1;
                t = (double)(i) * (2 * m + 1) / (double)(lim);
                s2 = ceil(s) - s;
                t1 = t - floor(t);

                is = (int)ceil(s) - 1;
                it = (int)floor(t) - 1;
                if (is < 0) is = 0;
                if (it < 0) it = 0;
                if (is >= n) is = n - 1;
                if (it >= n) it = n - 1;

                parz = 0;
                jstart = is + 1;
                jend_excl = it;
                if (jstart < 0) jstart = 0;
                if (jend_excl < 0) jend_excl = 0;
                if (jstart > n) jstart = n;
                if (jend_excl > n) jend_excl = n;
                if (jstart > jend_excl) jstart = jend_excl;
                for (int j = jstart; j < jend_excl; j++)
                    parz += y[j];

                a[i - 1] = parz + s2 * y[is] + t1 * y[it];
            }
            *dim_a = lim;
        }
        else
        {
            new_k = floor(k);
            extra = k - new_k;
            c = (2 * m + 1) / (2 * new_k + 1 + 2 * extra);

            a = (double *)calloc((size_t)(2 * (int)new_k + 3), sizeof(double));
            if (a == NULL)
                return 0;

            t = extra * c + 1;
            t1 = t - floor(t);

            if (k < 0)
            {
                free(a);
                return 0;
            }

            {
                int it0 = (int)floor(t) - 1;
                int jend0 = (int)floor(t);
                if (it0 < 0) it0 = 0;
                if (it0 >= n) it0 = n - 1;
                if (jend0 < 0) jend0 = 0;
                if (jend0 > n) jend0 = n;

                parz = 0;
                for (int i = 0; i < it0; i++)
                    parz += y[i];

                a[0] = parz + t1 * y[it0];
            }

            for (int i = 2; i < (int)(2 * new_k) + 3; i++)
            {
                int is, it, jstart, jend_excl;
                s = extra * c + (double)(i - 2) * c + 1;
                t = extra * c + (double)(i - 1) * c;
                s2 = ceil(s) - s;
                t1 = t - floor(t);

                is = (int)ceil(s) - 1;
                it = (int)floor(t) - 1;
                if (is < 0) is = 0;
                if (it < 0) it = 0;
                if (is >= n) is = n - 1;
                if (it >= n) it = n - 1;

                parz = 0;
                jstart = is + 1;
                jend_excl = it;
                if (jstart < 0) jstart = 0;
                if (jend_excl < 0) jend_excl = 0;
                if (jstart > n) jstart = n;
                if (jend_excl > n) jend_excl = n;
                if (jstart > jend_excl) jstart = jend_excl;
                for (int j = jstart; j < jend_excl; j++)
                    parz += y[j];

                a[i - 1] = parz + s2 * y[is] + t1 * y[it];
            }

            t2 = ceil(t) - t;

            {
                int is_tail = (int)ceil(t) - 1;
                int jstart_tail = (int)ceil(t) - 1;
                if (is_tail < 0) is_tail = 0;
                if (is_tail >= n) is_tail = n - 1;
                if (jstart_tail < 0) jstart_tail = 0;
                if (jstart_tail > n) jstart_tail = n;

                parz = 0;
                for (int i = jstart_tail + 1; i < n; i++)
                    parz += y[i];

                a[(int)(2 * new_k) + 2] = parz + t2 * y[is_tail];
            }
            *dim_a = 2 * (int)new_k + 3;
        }
    }
    else
    { 
        dx = 0.01;
        f_appo = (double *)calloc((size_t)n, sizeof(double));
        if (f_appo == NULL)
            return 0;
        f_appoB = 1;
        
        double *xappo = (double *)calloc((size_t)((int)m + 1), sizeof(double));
        if (xappo == NULL) {
            if (f_appoB) {
                free(f_appo);
                f_appoB = 0;
            }
            return 0;
        }
        for (int i = 0; i < ((int)m + 1); i++) xappo[i] = (double)i;

        for (int i = 0; i < n; i++) f_appo[i] = y[i] / dx;
        if (k == 0) {
            free(xappo);
            free(f_appo);
            return 0;
        }
        dy = m * dx / k;
        (void)dy;

        lim = (int)floor(k); 
        if (lim <= 0) {
            free(xappo);
            if (f_appoB) {
                free(f_appo);
                f_appoB = 0;
            }
            return 0;
        }

        double *xappo_interp = (double *)calloc((size_t)lim, sizeof(double));
        if (xappo_interp == NULL) {
            free(xappo);
            if (f_appoB) {
                free(f_appo);
                f_appoB = 0;
            }
            return 0;
        }
        for (int i = 0; i < lim; i++) xappo_interp[i] = (double)i * m / k;
        
        double *f_appo_interp = interp_linear(1, (int)m + 1, xappo, &f_appo[(int)m], lim, xappo_interp);
        if (f_appo_interp == NULL) {
            free(xappo_interp);
            free(xappo);
            if (f_appoB) {
                free(f_appo);
                f_appoB = 0;
            }
            return 0;
        }

        a = (double *)calloc((size_t)(2 * lim - 1), sizeof(double));
        if (a == NULL) {
            free(f_appo_interp);
            free(xappo_interp);
            free(xappo);
            if (f_appoB) {
                free(f_appo);
                f_appoB = 0;
            }
            return 0;
        }
        *dim_a = 2 * lim - 1; 

        for (int i = 0; i < lim - 1; i++) a[i] = f_appo_interp[lim - 1 - i];
        memcpy(a + lim - 1, f_appo_interp, (size_t)lim * sizeof(double));

        free(f_appo_interp);
        free(xappo_interp);
        free(xappo);

        Norm1 = 0;
        for (int i = 0; i < *dim_a; i++)
            Norm1 += fabs(a[i]);

        if (fabs(Norm1 - 1) > tol)
        {
            if (Norm1 == 0) {
                free(a);
                if (f_appoB)
                    free(f_appo);
                return 0;
            }

            for (int i = 0; i < *dim_a; i++)
                a[i] /= Norm1;
        }
    }

    if (f_appoB)
        free(f_appo);

    return a;
}

//********************** Maxmins_v3_8 *******************************
Maxmins Maxmins_v3_8(double *f, unsigned int N, double tol)
{
    unsigned int *Maxs;
    unsigned int *Mins;
    double *df;
    double *f1, *z = NULL;  // z never allocated here; initialize for safety
    unsigned int N_old, last_df = 0;
    double fappo1, fappo2, fappo3;
    int h, cmaxs, cmins, c;
    int posc = 0;
    int ctot, mac = 0, mic = 0;
    Maxmins valori;
    bool dfB = 0, f1B = 0;
    bool MMB = 0, MmB = 0;

    (void)z;
    (void)mac;
    (void)mic;

    valori.nout = 0;
    valori.maxmins = NULL;

    if (f == NULL || N < 2)
        return valori;

    Maxs = (unsigned int *)malloc(sizeof(unsigned int) * N);
    Mins = (unsigned int *)malloc(sizeof(unsigned int) * N);
    df   = (double *)malloc(sizeof(double) * (N - 1));

    if (Maxs == NULL || Mins == NULL || df == NULL)
    {
        if (Maxs) free(Maxs);
        if (Mins) free(Mins);
        if (df) free(df);
        return valori;
    }

    MMB = MmB = dfB = 1;

    for (int x = 0; x < N - 1; x++)
        df[x] = f[x + 1] - f[x];

    h = 1;
    while ((h < N) && (f[h - 1] != 0.0) && (fabs(df[h - 1] / f[h - 1]) <= tol))
        h++;

    if (dfB)
    {
        free(df);
        dfB = 0;
    }

    cmaxs = 0;
    cmins = 0;
    c = 0;
    N_old = N;

    if ((f1 = (double *)malloc(sizeof(double) * (N + h))) != NULL)
        f1B = 1;
    if ((df = (double *)malloc(sizeof(double) * (N + h))) != NULL)
        dfB = 1;
    if (!f1B || !dfB)
    {
        if (f1B) free(f1);
        if (dfB) free(df);
        if (MmB) free(Mins);
        if (MMB) free(Maxs);
        valori.nout = 0;
        valori.maxmins = NULL;
        return valori;
    }

    for (int i = 0; i < N; i++)
        f1[i] = f[i];
    for (int i = 0; i < h; i++)
        f1[N + i] = f[i + 1];

    N = N + h;

    for (int x = 0; x < N - 1; x++)
        df[x] = f1[x + 1] - f1[x];

    for (int i = h - 1; i < N - 2; i++)
    { 
        fappo2 = ((df[i] * df[i + 1]) / (fabs(f1[i]) * fabs(f1[i])));

        if ((fappo2 <= tol) && (fappo2 >= -tol))
        {
            fappo3 = df[i] / fabs(f1[i]);
            if (fappo3 < -tol)
            {
                last_df = -1;
                posc = i;
            }
            else if (fappo3 > tol)
            {
                last_df = +1;
                posc = i;
            }
            else if (df[i] == 0)
            {
                last_df = 0;
                posc = i;
            }

            c++;
            fappo1 = df[i + 1] / fabs(f1[i]);
            if (fappo1 < -tol)
            {
                if ((last_df == +1) || (last_df == 0))
                    Maxs[cmaxs++] = ((posc + (int)floor((double)(c - 1) / 2) + 1) % (N_old - 1));
                c = 0;
            }
            if (fappo1 > tol)
            {
                if ((last_df == -1) || (last_df == 0))
                    Mins[cmins++] = ((posc + (int)floor((double)(c - 1) / 2) + 1) % (N_old - 1));
                c = 0;
            }
        }

        if (fappo2 < -tol)
        {
            fappo3 = df[i] / fabs(f1[i]);
            fappo1 = df[i + 1] / fabs(f1[i]);

            if ((fappo3 < -tol) && (fappo1 > tol))
            {
                Mins[cmins++] = ((i + 1) % (N_old - 1));
                last_df = -1;
            }
            else if ((fappo3 > tol) && (fappo1 < -tol))
            {
                Maxs[cmaxs++] = ((i + 1) % (N_old - 1));
                last_df = +1;
            }
        }
    }

    if (c > 0)
    {
        if ((cmins > 0) && (Mins[cmins] == 0))
            Mins[cmins] = N;
        if ((cmaxs > 0) && (Maxs[cmaxs] == 0))
            Maxs[cmaxs] = N;
    }

    ctot = cmaxs + cmins;
    valori.nout = ctot;
    if (ctot > 0)
    {
        valori.maxmins = (unsigned int *)malloc(sizeof(unsigned int) * ctot);
        if (valori.maxmins == NULL)
        {
            valori.nout = 0;
            if (f1B) free(f1);
            if (dfB) free(df);
            if (MmB) free(Mins);
            if (MMB) free(Maxs);
            return valori;
        }

        for (int i = 0; i < cmaxs; i++) valori.maxmins[i] = Maxs[i];
        for (int i = 0; i < cmins; i++) valori.maxmins[i + cmaxs] = Mins[i];
        qsort(valori.maxmins, ctot, sizeof(unsigned int), comp);
    }

    if (f1B) free(f1);
    if (dfB) free(df);
    if (MmB) free(Mins);
    if (MMB) free(Maxs);

    // z is never allocated in this function; do not free it

    return valori;
}

//********************** FIF_v2_1 *******************************
static void free_maxmins_struct(Maxmins *mm)
{
    if (mm != NULL && mm->maxmins != NULL)
    {
        free(mm->maxmins);
        mm->maxmins = NULL;
        mm->nout = 0;
    }
}

Fif_t FIF_v2_1(double *f, int N, int *maxIMF)
{
    int numIMF = 0, completedIMF = 0, dim_a, nh;
    unsigned int k_pp;
    int in_step, ext_sig;
    int posF; 
    double df;
    int Nza, Nxs, N_old, N_r;
    int cont = 0;
    double rho, theta, thetan;
    double Norm1 = 1, m, logM = 0;
    double SD, norm_n, norm_d;
    double *h = NULL;
    double *mask = NULL, *mappo = NULL;
    double *ifftA = NULL, *fappo = NULL, *z = NULL;
    int n;

    n = *maxIMF;

    FILE *fpin = NULL;
    unsigned int conta = 0, size;
    bool hB, ifftAB, fappoB, maskB, mappoB;

    fif_complex *fftH = NULL, *ffth_new = NULL, *ffth_old = NULL;
    Fif_t IMF, *IMFappo, *next;
    double tol = 10e-12;

    (void)nh;
    (void)cont;
    (void)rho;
    (void)theta;
    (void)thetan;
    (void)z;
    (void)fpin;
    (void)conta;

    memset(&IMF, 0, sizeof(IMF));

    if (maxIMF == NULL)
        return IMF;

    hB = ifftAB = fappoB = maskB = mappoB = 0;
    free_maxmins_struct(&MaxMin);

    IMFappo = &IMF;
    size = sizeM;

    Norm1 = 0;
    for (int i = 0; i < N; i++)
        Norm1 += fabs(f[i]);

    if (Norm1 == 0)
        goto cleanup_success;

    for (int i = 0; i < N; i++)
        f[i] /= Norm1;

    N_r = 0;
    fappo = (double *)calloc((size_t)(N + 11), sizeof(double));
    if (fappo == NULL)
        goto cleanup_success;
    fappoB = 1;

    for (int i = 0; i < N; i++)
    {
        if (fabs(f[i]) > tol)
            fappo[N_r++] = f[i];
    }
    if (N_r == 0)
        goto cleanup_success;

    for (int i = 0; i < 10; i++)
        fappo[N_r + i] = fappo[i];

    free_maxmins_struct(&MaxMin);
    MaxMin = Maxmins_v3_8(fappo, N_r + 10, tol);

    if (fappoB)
    {
        free(fappo);
        fappo = NULL;
        fappoB = 0;
    }

    k_pp = MaxMin.nout;
    if (k_pp == 0)
        goto cleanup_success;

/*    for (int i = 0; i < MaxMin.nout; i++)
        printf("%d ", MaxMin.maxmins[i]);
    printf("\n");  Spammy printout removed  CC 2/21/26 */

    numIMF = 0;
    cont = 0;

    while ((completedIMF < n) && (k_pp > ExtPoints))
    {
        numIMF = completedIMF + 1;
        cont++;
        SD = 1;

        h = (double *)malloc(sizeof(double) * N);
        if (h == NULL)
            goto cleanup_success;
        hB = 1;
        memcpy(h, f, sizeof(double) * N);
        nh = N;

        m = round(2 * N_r / (k_pp + 1) * Xi);

        if (numIMF > 1)
        {
            if (m <= logM)
            {
                m = ceil(logM * 1.1);
            }
        }

        logM = m;
        IMFappo->stats.logM = m;

        if (maskB)
        {
            maskB = 0;
            free(mask);
            mask = NULL;
        }

        mask = getMask(MM, size, m, &dim_a, tol);
        if (mask == NULL)
            goto cleanup_success;
        maskB = 1;

        ext_sig = 0;
        N_r = N;

        if (N < dim_a)
        {
            ext_sig = 1;
            Nxs = (int)ceil((float)dim_a / (float)N);
            N_old = N;
            if ((Nxs % 2) == 0)
                Nxs++;
            N_r = N * Nxs;

            if (fappoB)
            {
                free(fappo);
                fappo = NULL;
                fappoB = 0;
            }

            fappo = (double *)malloc(sizeof(double) * N_r);
            if (fappo == NULL)
                goto cleanup_success;
            fappoB = 1;

            for (int i = 0; i < Nxs; i++)
                memcpy(fappo + N * i, h, N * sizeof(double));

            if (hB)
            {
                free(h);
                h = NULL;
                hB = 0;
            }

            h = (double *)malloc(sizeof(double) * N_r);
            if (h == NULL)
                goto cleanup_success;
            hB = 1;
            memcpy(h, fappo, sizeof(double) * N_r);

            if (fappoB)
            {
                free(fappo);
                fappo = NULL;
                fappoB = 0;
            }
        }

        Nza = N_r - dim_a;

        if (mappoB)
        {
            free(mappo);
            mappo = NULL;
            mappoB = 0;
        }

        mappo = (double *)calloc((size_t)N_r, sizeof(double));
        if (mappo == NULL)
            goto cleanup_success;
        mappoB = 1;

        if ((Nza % 2) == 0)
        {
            memcpy(mappo, mask + dim_a / 2, (dim_a / 2) * sizeof(double));
            memcpy(mappo + (N_r - dim_a / 2), mask, (dim_a / 2) * sizeof(double));
        }
        else
        {
            memcpy(mappo, mask + dim_a / 2, (dim_a / 2 + 1) * sizeof(double));
            memcpy(mappo + (N_r - dim_a / 2), mask, (dim_a / 2) * sizeof(double));
        }

        if (ifftAB)
        {
            free(ifftA);
            ifftA = NULL;
            ifftAB = 0;
        }

        fflush(NULL);
        ifftA = realFFT(mappo, N_r);
        fflush(NULL);
        if (ifftA == NULL)
            goto cleanup_success;
        ifftAB = 1;

        if (mappoB)
        {
            free(mappo);
            mappo = NULL;
            mappoB = 0;
        }

        ffth_new = (fif_complex *)malloc(sizeof(fif_complex) * N_r);
        ffth_old = (fif_complex *)malloc(sizeof(fif_complex) * N_r);
        fftH = fft_dir(h, N_r);
        if (ffth_new == NULL || ffth_old == NULL || fftH == NULL)
            goto cleanup_success;

        posF = 0;
        df = -1.;
        while ((posF < N_r - 1) && (df < 0))
        {
            df = ifftA[posF + 1] - ifftA[posF];
            posF++;
        }
        posF--;

        IMFappo->stats.posF = posF;
        IMFappo->stats.valF = ifftA[posF];

        for (int i = 0; i < N_r; i++)
        {
            ifftA[i] -= IMFappo->stats.valF;
            if (ifftA[i] < 0)
                ifftA[i] = 0;
        }

        in_step = 0;
        while ((SD > DELTA) && (in_step < MaxInner))
        {
            in_step += 1;

            norm_n = norm_d = 0;
            for (int i = 0; i < N_r; i++)
            {
                ffth_new[i].re = pow((1 - ifftA[i]), in_step) * fftH[i].re;
                ffth_new[i].im = pow((1 - ifftA[i]), in_step) * fftH[i].im;

                ffth_old[i].re = pow((1 - ifftA[i]), in_step - 1) * fftH[i].re;
                ffth_old[i].im = pow((1 - ifftA[i]), in_step - 1) * fftH[i].im;

                norm_n += (ffth_new[i].re - ffth_old[i].re) * (ffth_new[i].re - ffth_old[i].re) +
                          (ffth_new[i].im - ffth_old[i].im) * (ffth_new[i].im - ffth_old[i].im);
                norm_d += ffth_old[i].re * ffth_old[i].re + ffth_old[i].im * ffth_old[i].im;
            }

            SD = norm_n / norm_d;
        }

        if (hB)
        {
            free(h);
            h = NULL;
            hB = 0;
        }

        h = fft_inv(ffth_new, N_r);

        free(fftH);
        free(ffth_new);
        free(ffth_old);
        fftH = ffth_new = ffth_old = NULL;

        if (h == NULL)
            goto cleanup_success;
        hB = 1;

        if (ext_sig)
        {
            N = N_old;
            int j = 0;
            for (int i = (N * (Nxs - 1) / 2); i < (N * ((Nxs - 1) / 2 + 1)); i++)
                h[j++] = h[i];

            if (j == N)
            {
                double *tmp = (double *)realloc(h, sizeof(double) * N);
                if (tmp != NULL)
                    h = tmp;
            }
        }

        IMFappo->stats.in_step = in_step;
        IMFappo->dati = (double *)malloc(sizeof(double) * N);
        if (IMFappo->dati == NULL)
            goto cleanup_success;

        for (int i = 0; i < N; i++)
            IMFappo->dati[i] = h[i] * Norm1;

        for (int i = 0; i < N; i++)
            f[i] -= h[i];

        next = (Fif_t *)malloc(sizeof(Fif_t));
        if (next == NULL)
        {
            free(IMFappo->dati);
            IMFappo->dati = NULL;
            goto cleanup_success;
        }
        memset(next, 0, sizeof(Fif_t));
        IMFappo->next = next;           
        IMFappo = next;
        completedIMF++;

        N_r = 0;

        fappo = (double *)calloc((size_t)(N + 10), sizeof(double));
        if (fappo == NULL)
            goto cleanup_success;
        fappoB = 1;

        for (int i = 0; i < N; i++)
        {
            if (fabs(f[i]) > tol)
                fappo[N_r++] = f[i];
        }

        if (N_r == 0)
            goto cleanup_success;

        memcpy(fappo + N_r, fappo, 10 * sizeof(double));

        free_maxmins_struct(&MaxMin);
        if (logM >= 20)
        {
            double *sappo = (double *)malloc(sizeof(double) * (N + 10));
            if (sappo == NULL)
                goto cleanup_success;
            sappo[0] = fappo[0]; 
            sappo[1] = (fappo[0] + fappo[1]) / 2;
            sappo[2] = (fappo[0] + fappo[1] + fappo[2] + fappo[3]) / 4;
            sappo[3] = (fappo[0] + fappo[1] + fappo[2] + fappo[3] + fappo[4] + fappo[5]) / 6;
            sappo[4] = (fappo[0] + fappo[1] + fappo[2] + fappo[3] + fappo[4] +
                        fappo[5] + fappo[6] + fappo[7]) / 8;

            sappo[N_r + 9] = fappo[N_r + 9]; 
            sappo[N_r + 8] = (fappo[N_r + 9] + fappo[N_r + 8]) / 2;
            sappo[N_r + 7] = (fappo[N_r + 9] + fappo[N_r + 8] + fappo[N_r + 7] + fappo[N_r + 6]) / 4;
            sappo[N_r + 6] = (fappo[N_r + 9] + fappo[N_r + 8] + fappo[N_r + 7] + fappo[N_r + 6] +
                              fappo[N_r + 5] + fappo[N_r + 4]) / 6;
            sappo[N_r + 5] = (fappo[N_r + 9] + fappo[N_r + 8] + fappo[N_r + 7] + fappo[N_r + 6] +
                              fappo[N_r + 5] + fappo[N_r + 4] + fappo[N_r + 3] + fappo[N_r + 2]) / 8;

            for (int i = 5; i < N_r + 5; i++)
            {
                sappo[i] = (fappo[i - 5] + fappo[i - 4] + fappo[i - 3] + fappo[i - 2] + fappo[i - 1] +
                            fappo[i] + fappo[i + 1] + fappo[i + 2] + fappo[i + 3] + fappo[i + 4]) / 10;
            }

            memcpy(fappo, sappo, sizeof(double) * (N_r + 10));
            MaxMin = Maxmins_v3_8(fappo, N_r + 10, tol);
            free(sappo);
        }
        else
        {
            MaxMin = Maxmins_v3_8(fappo, N_r + 10, tol);
        }

        k_pp = MaxMin.nout;

        if (fappoB)
        {
            free(fappo);
            fappo = NULL;
            fappoB = 0;
        }

        if (hB)
        {
            free(h);
            h = NULL;
            hB = 0;
        }

        if (k_pp == 0)
        {
            goto cleanup_success;
        }
    }

    // save residual in the current tail node
    IMFappo->dati = (double *)malloc(sizeof(double) * N);
    if (IMFappo->dati != NULL)
    {
        for (int i = 0; i < N; i++)
            IMFappo->dati[i] = f[i] * Norm1;
    }
    IMFappo->next = NULL;

cleanup_success:
    if (fftH != NULL) free(fftH);
    if (ffth_new != NULL) free(ffth_new);
    if (ffth_old != NULL) free(ffth_old);
    if (hB && h != NULL) free(h);
    if (maskB && mask != NULL) free(mask);
    if (mappoB && mappo != NULL) free(mappo);
    if (ifftAB && ifftA != NULL) free(ifftA);
    if (fappoB && fappo != NULL) free(fappo);
    free_maxmins_struct(&MaxMin);

    *maxIMF = completedIMF;
    fflush(NULL);
    return IMF;
}