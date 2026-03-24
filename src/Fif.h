#ifndef FIF_FIF_H

typedef struct Maxmins Maxmins;
struct Maxmins {
    int nout;
    unsigned int *maxmins;
};

#define FIF_FIF_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct statistiche statistiche;
typedef struct Fif_t Fif_t;

struct statistiche {
    double logM;
    int posF;
    double valF;
    int in_step;
};

struct Fif_t {
   statistiche stats;
   Fif_t * next;
   double * dati;
};

int estimatedIMF(double *h,int N,int *k_pp);
Fif_t FIF_v2_1(double *f, int N, int *n);

#ifdef __cplusplus
}
#endif

#endif
