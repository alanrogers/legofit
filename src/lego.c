/**
 * @file lego.c
 * @brief Simulate branch lengths
 */

#include "gptree.h"
#include "binary.h"
#include "jobqueue.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct TaskArg TaskArg;

// Number of site patterns
#define N_SITE_PAT 6

/** Data structure used by each thread */
struct TaskArg {
    /* rates of gene flow into E from N and from D */
    long double mN, mD;

    /* time in units of 2N0 generations */
    long double alpha, delta, zeta, kappa, lambda;
    long double lambda_max;

    // mutation rate per 2N0 generations
    double mutRate;

    /* population sizes relative to N0 */
    long double K_D, K_ND, K_N, K_XY;

    unsigned rng_seed;
    unsigned long nreps;

    /* Returned values */
    long unsigned sitePatCount[N_SITE_PAT];
    long double Ixnd; /* fraction of sites matching */
    long double oxnd; /* simulation averages */
    long double Ixn_; /* fraction of sites matching */
    long double oxn_; /* simulation averages */
    long double Ixn; /* fraction of sites matching */
    long double oxn; /* simulation averages */
    long double Ixy; /* fraction of sites matching */
    long double oxy; /* simulation averages */
    long double Iyn_; /* fraction of sites matching */
    long double oyn_; /* simulation averages */
    long double Iyn; /* fraction of sites matching */
    long double oyn; /* simulation averages */
    long double Iyd; /* fraction of sites matching */
    long double oyd; /* simulation averages */
    long double Iynd; /* fraction of sites matching */
    long double oynd; /* simulation averages */
    long double Iyd_;
    long double oyd_;
};

/* Prototypes */
void TaskArg_free(TaskArg *a);

int taskfun(void *varg);

TaskArg *TaskArg_new(const TaskArg *template, unsigned rng_seed);

void TaskArg_printParams(TaskArg *a, FILE *fp);

void doSim(unsigned long nreps, PopNode *rootPop, PopNode *A,
             PopNode *B, PopNode *C, PopNode *D, tipId_t tipId[], gsl_rng *rng,
             TaskArg *a);

double *XY_pat(double lambda, double kappa, double zeta,
           double K_XY, double K_ND, double m_XY);

double *N_pat(double lambda, double kappa, double delta,
          double K_N, double K_ND, double m_N);

double *D_pat(double lambda, double kappa, double alpha,
          double K_D, double K_ND, double m_D);

long double max(long double x, long double y);

double sum(double *x, int n);

long double max(long double x, long double y) {
    if(x > y){
        return x;
    } else {
        return y;
    }
}

double sum(double *x, int n) {
    int i;
    double tot = 0;
    for(i=0; i < n; i++){
        tot = tot + x[i];
    }
    return tot;
}

void TaskArg_printParams(TaskArg *a, FILE *fp) {
    fprintf(fp, "%-20s: %Lf\n", "alpha", a->alpha);
    fprintf(fp, "%-20s: %Lf\n", "delta", a->delta);
    fprintf(fp, "%-20s: %Lf\n", "zeta", a->zeta);
    fprintf(fp, "%-20s: %Lf\n", "kappa", a->kappa);
    fprintf(fp, "%-20s: %Lf\n", "lambda", a->lambda);
    fprintf(fp, "%-20s: %Lf\n", "K_D", a->K_D);
    fprintf(fp, "%-20s: %Lf\n", "K_N", a->K_N);
    fprintf(fp, "%-20s: %Lf\n", "K_ND", a->K_ND);
    fprintf(fp, "%-20s: %Lf\n", "K_XY", a->K_XY);
    fprintf(fp, "%-20s: %lu\n", "nreps", a->nreps);

    putc('\n', fp);
}

void doSim(unsigned long nreps, PopNode *rootPop, PopNode *A,
             PopNode *B, PopNode *C, PopNode *D, tipId_t tipId[], gsl_rng *rng,
             TaskArg *a) {
    unsigned long rep;
    unsigned i;
    Gene *root;

    /*
     * To be scrupulous, I should ignore iterations in which the
     * genealogy has a human outgroup--ones in which one human lineage
     * is closer to chimp than to other humans. However, these arise
     * only when the two human lineages survive without coalescing
     * into the Chuman population, and in that case the distribution
     * of branch lengths is the same no matter which topology we
     * get. So I'm just ignoring topology.
     */

    a->oxnd = 0;
    a->oxn = 0;
    a->oxn_ = 0;
    a->oxy = 0;
    a->oyn = 0;
    a->oyd = 0;
    a->oynd = 0;
    a->oyn_ = 0;

    /**
     * tipId array is in the following order:
     *      X: 1; Y: 2; N: 4; D: 8;
     *   long double oxnd = 1|4|8
     *   long double oxn = 1|4
     *   long double oxy = 1|2
     *   long double oyn = 2|4
     *   long double oyd = 2|8
     *   long double oynd = 2|4|8
     *
     * tipId_t tipIds[6] = {1|4|8, 1|4, 1|2, 2|4, 2|8, 2|4|8};
     */

    for(rep = 0; rep < nreps; ++rep) {
        /* put samples into populations */
        PopNode_clear(rootPop);
        PopNode_newGene(A, 0);
        PopNode_newGene(B, 1);
        PopNode_newGene(C, 2);
	    PopNode_newGene(D, 3);

        /* coalescent simulation */
        root = PopNode_coalesce(rootPop, rng);
        assert(root);
	    /* identifies the tipId in the function, and assigns from there */
        double bl[N_SITE_PAT]; // branch lengths
        double blTot = 0.0;
        for(i=0; i < N_SITE_PAT; ++i)
            blTot += bl[i] = Gene_getRightLen(root, tipId[i]);

        // Probability that a mutation occurs on this tree
        double prMut = 1.0 - exp(-blTot*a->mutRate);

        double u = gsl_rng_uniform(rng);
        if( u < prMut) {
            // A mutation happened. Figure out where.
            double cum[N_SITE_PAT];
            cum[0] = bl[0];

            // Cumulative branch lengths
            for(i=1; i< N_SITE_PAT; ++i)
                cum[i] = bl[i] + cum[i-1];

            u *= cum[N_SITE_PAT-1];
            i = Dbl_first_geq(u, N_SITE_PAT, cum);
            assert(i < N_SITE_PAT);
            a->sitePatCount[i] += 1;
        }

        a->oxnd += bl[0];
        a->oxn += bl[1];
        a->oxy += bl[2];
        a->oyn += bl[3];
        a->oyd += bl[4];
        a->oynd += bl[5];

        Gene_free(root);
    }
    a->oxn_ += a->oxnd + a->oxn;
    a->oyn_ += a->oyn + a->oynd;

    a->oxnd /= nreps;
    a->oxn /= nreps;
    a->oxy /= nreps;
    a->oyn /= nreps;
    a->oyd /= nreps;
    a->oynd /= nreps;
    a->oxn_ /= nreps;
    a->oyn_ /= nreps;
}

double *XY_pat(double lambda, double kappa, double zeta,
           double K_XY, double K_ND, double m_XY) {
    static double ret[7];

    double p_xn;
    double p_xy;
    double p_yn_;
    double p_yn;
    double p_yd_;
    double p_yd;
    double p_ynd;

    double z = zeta;
    double l = lambda;

    /**
     * Calculate all of the terms necessary for
     * fully calculating the probability and
     * expected branch length vectors.
     */

     // Just survival functions for these intervals
     double x_kl[1] = {exp(-(lambda - kappa) / K_ND)};
     double x_zl[1] = {exp(-(lambda - zeta) / K_XY)};

     // Expected branch length vectors for XY and ND
     double m_zl[1] = {K_XY * (1. - x_zl[0])};

     /**
      * 3 cases in Anc:
      * 4LOD = x, y --> Anc :: n, d --> Anc
      * 3LOD = xy --> Anc :: n, d --> Anc
      *      = x, y --> Anc :: nd --> Anc
      * 2LOD = xy --> Anc :: nd --> Anc
      */

    // Branch length vectors associated with prior initial vectors
    double m_2l[1] = {1.};
    double m_3l[2] = {1., 1./3.};
    double m_4l[3] = {1., 1./3., 1./6.};

    p_xn = (m_XY * x_zl[0] * (m_3l[0] / 3.));
    p_xy = (
        m_XY * (
            (l - z - m_zl[0]) +
            x_zl[0] * m_3l[0] / 3. +
            (1. - x_zl[0]) * m_2l[0]
        )
    );
    p_yn_ = (
        m_XY * (
            x_zl[0] * (
                m_3l[0] / 3.)
            )
        );
    p_yn = (m_XY * (x_zl[0] * x_kl[0] * (1./6. * (m_4l[1] + m_4l[0] / 3.)
                                         + 1./6. * 1./3. * m_4l[0])));
    p_yd_ = (m_XY * (x_zl[0] * (m_3l[0] / 3.)));
    p_yd = (m_XY * (x_zl[0] * x_kl[0] * (1./6. * (m_4l[1] + m_4l[0] / 3.)
                                         + 1./6. * 1./3. * m_4l[0])));
    p_ynd = (m_XY * (x_zl[0] * x_kl[0] * (1./2. * 1./3. * m_4l[0])
                     + x_zl[0] * (1 - x_kl[0]) * (m_3l[0] / 3.)));

    ret[0] = p_xn;
    ret[1] = p_xy;
    ret[2] = p_yn;
    ret[3] = p_yd;
    ret[4] = p_ynd;
    ret[5] = p_yn_;
    ret[6] = p_yd_;

    return ret;
}

double *N_pat(double lambda, double kappa, double delta,
          double K_N, double K_ND, double m_N) {
    static double ret[7];

    double p_xn;
    double p_xy;
    double p_yn_;
    double p_yn;
    double p_yd_;
    double p_yd;
    double p_ynd;

    double d = delta;
    double k = kappa;
    double l = lambda;

    /**
     * Calculate all of the terms necessary for
     * fully calculating the probability and
     * expected branch length vectors.
     */
     // Probability vectors for N
     double x_dk[1] = {exp(-(kappa - delta) / K_N)};

     // Expected branch length for N
     double m_dk[1] = {K_N * (1. - exp(-(kappa - delta) / K_N))};

     double x_2kl[1] = {exp(-(lambda - kappa) / K_ND)};
     double x_3kl[2] = {
        3./2. * exp(-(lambda - kappa) / K_ND)
        - 3./2. * exp(-3. * (lambda - kappa) / K_ND),
        exp(-3. * (lambda - kappa) / K_ND)};

     // Expected branch length vectors for ND
     double m_2kl[1] = {K_ND * (1. - exp(-(lambda - kappa) / K_ND))};
     double m_3kl[2] = {
        K_ND * (1./2. * exp(-3. * (lambda - kappa) / K_ND)
                - 3./2. * exp(-(lambda - kappa) / K_ND) + 1.),
        K_ND * 1./3. * (1 - exp(-3. * (lambda-kappa) / K_ND))};

     /**
      * 3 cases in Anc:
      * 4LOD = x --> Anc :: y, n, d --> Anc
      * 3LOD = x --> Anc :: yn, d or y, nd or yd, n --> Anc
      * 2LOD = x --> Anc :: ynd --> Anc
      */

    // Branch length vectors associated with prior initial vectors
    double m_2l[1] = {1.};
    double m_3l[2] = {1., 1./3.};
    double m_4l[3] = {1., 1./3., 1./6.};

    p_xn = (
        m_N * (
            x_dk[0] * x_2kl[0] * (m_3l[0] / 3.)
        )
    );

    p_xy = p_xn;

    p_yn_ = (
        m_N * (
            (k - d - m_dk[0]) +
            x_dk[0] * (
                (l - k - m_2kl[0]) +
                x_2kl[0] * (
                    m_3l[0] / 3.) +
                (1 - x_2kl[0]) * (
                    m_2l[0])
                ) +
            (1 - x_dk[0]) * (
                (l - k) +
                m_2l[0])
        )
    );

    p_yn = (
        m_N * (
            (k - d - m_dk[0]) +
            x_dk[0] * (
                m_3kl[0] / 3. +
                x_3kl[1] * (
                    1./6. * (m_4l[1] + 1./3. * m_4l[0]) +
                    1./6. * 1./3. * m_4l[0]) +
                x_3kl[0] * (
                    1./3. * (m_3l[1]  + 1./3. * m_3l[0]))
            ) +
            (1 - x_dk[0]) * (
                m_2kl[0] +
                x_2kl[0] * (
                    m_3l[1] + m_3l[0] / 3.)
            )
        )
    );

    p_yd_ = (
        m_N * (
            (l - k - m_2kl[0]) +
            x_2kl[0] * (m_3l[0] / 3.) +
            (1 - x_2kl[0]) * m_2l[0]
        )
    );

    p_yd = (
        m_N * (
            x_dk[0] * (
                m_3kl[0] / 3. +
                x_3kl[1] * (
                    1./6. * (m_4l[1] + 1./3. * m_4l[0]) +
                    1./6. * 1./3. * m_4l[0]) +
                x_3kl[0] * (
                    1./3. * (m_3l[1] + 1./3. * m_3l[0]))
            )
        )
    );

    p_ynd = (
        m_N * (
            (1. - x_dk[0]) * (
                (l - k - m_2kl[0]) +
                x_2kl[0] * (m_3l[0] / 3.) +
                (1. - x_2kl[0]) * m_2l[0]
                ) +
            x_dk[0] * (
                (l - k - sum(m_3kl, 2)) +
                x_3kl[1] * (1./2. * 1./3. * m_4l[0]) +
                x_3kl[0] * (1./3. * m_3l[0]) +
                (1 - sum(x_3kl, 2)) * m_2l[0]
            )
        )
    );

    ret[0] = p_xn;
    ret[1] = p_xy;
    ret[2] = p_yn;
    ret[3] = p_yd;
    ret[4] = p_ynd;
    ret[5] = p_yn_;
    ret[6] = p_yd_;

    return ret;
}

double *D_pat(double lambda, double kappa, double alpha,
          double K_D, double K_ND, double m_D) {
    static double ret[7];

    double p_xn;
    double p_xy;
    double p_yn_;
    double p_yn;
    double p_yd_;
    double p_yd;
    double p_ynd;

    double a = alpha;
    double k = kappa;
    double l = lambda;

    /**
     * Calculate all of the terms necessary for
     * fully calculating the probability and
     * expected branch length vectors.
     */
     // Probability vectors for D
     double x_ak[1] = {exp(-(kappa - alpha) / K_D)};

     // Expected branch length for D
     double m_ak[1] = {K_D * (1. - exp(-(kappa - alpha) / K_D))};

     double x_2kl[1] = {exp(-(lambda - kappa) / K_ND)};
     double x_3kl[2] = {
        3./2. * exp(-(lambda - kappa) / K_ND)
        - 3./2. * exp(-3. * (lambda - kappa) / K_ND),
        exp(-3. * (lambda - kappa) / K_ND)};

     // Expected branch length vectors for ND
     double m_2kl[1] = {K_ND * (1. - exp(-(lambda - kappa) / K_ND))};
     double m_3kl[2] = {
        K_ND * (1./2. * exp(-3. * (lambda - kappa) / K_ND)
                - 3./2. * exp(-(lambda - kappa) / K_ND) + 1.),
        K_ND * 1./3. * (1 - exp(-3. * (lambda-kappa) / K_ND))};

     /**
      * 3 cases in Anc:
      * 4LOD = x --> Anc :: y, n, d --> Anc
      * 3LOD = x --> Anc :: yn, d or y, nd or yd, n --> Anc
      * 2LOD = x --> Anc :: ynd --> Anc
      */
    // Branch length vectors associated with prior initial vectors
    double m_2l[1] = {1.};
    double m_3l[2] = {1., 1./3.};
    double m_4l[3] = {1., 1./3., 1./6.};

    p_xn = (
        m_D * (
            x_2kl[0] * (m_3l[0] / 3.)
        )
    );

    p_xy = p_xn;

    p_yn_ = (
        m_D * (
            (l - k - m_2kl[0]) +
            x_2kl[0] * (m_3l[0] / 3.) +
            (1 - x_2kl[0]) * m_2l[0]
        )
    );

    p_yn = (
        m_D * (
            x_ak[0] * (
                m_3kl[0] / 3. +
                x_3kl[1] * (
                    1./6. * (m_4l[1] + 1./3. * m_4l[0]) +
                    1./6. * 1./3. * m_4l[0]) +
                x_3kl[0] * (
                    1./3. * (m_3l[1] + 1./3. * m_3l[0]))
            )
        )
    );

    p_yd_ = (
        m_D * (
            (k - a - m_ak[0]) +
            x_ak[0] * (
                (l - k - m_2kl[0]) +
                x_2kl[0] * (
                    m_3l[0] / 3.) +
                (1 - x_2kl[0]) * m_2l[0]) +
            (1 - x_ak[0]) * (
                (l - k) +
                m_2l[0])
        )
    );

    p_yd = (
        m_D * (
            (k - a - m_ak[0]) +
            x_ak[0] * (
                m_3kl[0] / 3. +
                x_3kl[1] * (
                    1./6. * (m_4l[1] + 1./3. * m_4l[0]) +
                    1./6. * 1./3. * m_4l[0]) +
                x_3kl[0] * (
                    1./3. * (m_3l[1] + 1./3. * m_3l[0]))
            ) +
            (1 - x_ak[0]) * (
                m_2kl[0] +
                x_2kl[0] * (
                    m_3l[1] + m_3l[0] / 3.)
            )
        )
    );

    p_ynd = (
        m_D * (
            (1. - x_ak[0]) * (
                (l - k - m_2kl[0]) +
                x_2kl[0] * (m_3l[0] / 3.) +
                (1. - x_2kl[0]) * m_2l[0]
                ) +
            x_ak[0] * (
                (l - k - sum(m_3kl, 2)) +
                x_3kl[1] * (1./2. * 1./3. * m_4l[0]) +
                x_3kl[0] * (1./3. * m_3l[0]) +
                (1 - sum(x_3kl, 2)) * m_2l[0]
            )
        )
    );

    ret[0] = p_xn;
    ret[1] = p_xy;
    ret[2] = p_yn;
    ret[3] = p_yd;
    ret[4] = p_ynd;
    ret[5] = p_yn_;
    ret[6] = p_yd_;

    return ret;
}

/**
 * Construct a new TaskArg by copying a template, but then assign
 * a distinct random number seed.
 */
TaskArg *TaskArg_new(const TaskArg *template, unsigned rng_seed) {
    TaskArg *a = malloc(sizeof(TaskArg));
    checkmem(a, __FILE__, __LINE__);

    memcpy(a, template, sizeof(TaskArg));
    a->rng_seed = rng_seed;

    memset(a->sitePatCount, 0, sizeof(a->sitePatCount));

    return a;
}

/** TaskArg destructor */
void TaskArg_free(TaskArg *a) {
    free(a);
}

/** function run by each thread */
int taskfun(void *varg) {
    TaskArg *a = (TaskArg *) varg;
    /**
     * Y0 is the tip, Y1 is what remains after alpha, and Y2 after delta.
     */

    PopNode *X, *Y0, *Y1, *Y2, *N, *D, *ND,
        *XY, *AEND;

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, a->rng_seed);

    /* perturb parameters */
#if 1
    if(a->mD > 0.0)
        a->mD = perturb_ratio(a->mD, rng);

    a->K_D = perturb_ratio(a->K_D, rng);
    a->K_N = perturb_ratio(a->K_N, rng);
    a->K_ND = perturb_ratio(a->K_ND, rng);
    a->K_XY = perturb_ratio(a->K_XY, rng);
#endif

#if 1
    /* zeta is guaranteed to be less than kappa in reality so
     * it should be fine to constrain it this way in simulation */
    a->alpha = perturb_interval(a->alpha, 0.5*a->alpha, a->zeta, rng);
    a->delta = perturb_interval(a->delta, 0.5*a->delta, a->zeta, rng);
    /* because alpha and delta are no longer fixed in order
     * we have to constrain zeta so that it is larger than whichever
     * is largest */
    a->zeta = perturb_interval(a->zeta, max(a->alpha, a->delta), a->kappa, rng);
    a->kappa = perturb_interval(a->kappa, a->zeta, a->lambda, rng);
    a->lambda = perturb_interval(a->lambda, a->kappa, a->lambda_max, rng);
#endif
    assert(a->alpha < a->kappa);
    assert(a->delta < a->kappa);
    assert(a->zeta < a->lambda);
    assert(a->kappa < a->lambda);

    X = PopNode_new(1.0, 0.0, a->zeta);
    Y0 = PopNode_new(1.0, 0.0, a->alpha);
    /* Tree should fit whether alpha < delta or
     * delta < alpha. This should ensure that is
     * true.
     *
     * I'm not sure what to do if they are equal,
     * but in this case it will just default to
     * the else statement. */
    if(a->alpha < a->delta){
        Y1 = PopNode_new(1.0, a->alpha, a->delta);
        Y2 = PopNode_new(1.0, a->delta, a->zeta);
    } else {
        Y1 = PopNode_new(1.0, a->delta, a->alpha);
        Y2 = PopNode_new(1.0, a->alpha, a->zeta);
    }
    N = PopNode_new(a->K_N, a->delta, a->kappa);
    D = PopNode_new(a->K_D, a->alpha, a->kappa);
    ND = PopNode_new(a->K_ND, a->kappa, a->lambda);
    XY = PopNode_new(a->K_XY, a->zeta, a->lambda);
    AEND = PopNode_new(1.0, a->lambda, DBL_MAX);

    if(a->alpha < a->delta) {
        /* Y0 is a mixture: mD * D + (1 - mD) * Y1 */
        PopNode_mix(Y0, a->mD, D, Y1);
        /* Y1 is a mixture: mN * N + (1 - mN) * Y2 */
        PopNode_mix(Y1, a->mN, N, Y2);
    } else {
        /* Y0 is a mixture: mN * N + (1 - mN) * Y1 */
        PopNode_mix(Y0, a->mN, N, Y1);
        /* Y1 is a mixture: mD * D + (1 - mD) * Y2 */
        PopNode_mix(Y1, a->mD, D, Y2);
    }

    PopNode_join(XY, X, Y2);
    PopNode_join(ND, N, D);
    PopNode_join(AEND, XY, ND);

    /**
     * This section generates nreps trees, with root population, and 1..3.
     * tipId array is in the following order:
     *      X: 1; Y: 2; N: 4; D: 8;
     *   long double oxnd = 1|4|8
     *   long double oxn = 1|4
     *   long double oxy = 1|2
     *   long double oyn = 2|4
     *   long double oyd = 2|8
     *   long double oynd = 2|4|8
     */
    tipId_t tipIds[6] = {1|4|8, 1|4, 1|2, 2|4, 2|8, 2|4|8};
    doSim(a->nreps, AEND, X, Y0, N, D, tipIds, rng, a);

    /* Sample paths */
    long double m_XY = (1.0 - a->mN - a->mD);
    long double m_N = a->mN;
    long double m_D = a->mD;

    /* Full estimate */
    /**
     * We'll be constructing this in sections, with XY, N, and D
     * having separate contributions.
     */
    double *xy_pat = XY_pat(a->lambda, a->kappa, a->zeta, a->K_XY, a->K_ND, m_XY);
    double *n_pat = N_pat(a->lambda, a->kappa, a->delta, a->K_N, a->K_ND, m_N);
    double *d_pat = D_pat(a->lambda, a->kappa, a->alpha, a->K_D, a->K_ND, m_D);

    /*
    double p_xn = (
        a->mN * exp(-(a->kappa - a->delta) / a->K_N) *
        exp(-(a->lambda - a->kappa) / a->K_ND) +
        a->mD * exp(-(a->lambda - a->kappa) / a->K_ND) +
        (1. - a->mN - a->mD) * exp(-(a->lambda - a->zeta) / a->K_XY)) / 3.;

    double p_yn_ = p_xn + (a->mN * (
        a->lambda - a->delta + (1. - a->K_N) * (1. - exp(
            -(a->kappa - a->delta) / a->K_N)) +
        (1. - a->K_ND) * exp(-(a->kappa - a->delta) / a->K_N) *
        (1. - exp(-(a->lambda - a->kappa) / a->K_ND)) +
        a->mD * (a->lambda - a->kappa + (1. - a->K_ND) * (
            1. - exp(-(a->lambda - a->kappa) / a->K_ND)))));
double p_yd_ = p_xn + (a->mD * (
        a->lambda - a->alpha + (1. - a->K_N) * (1. - exp(
            -(a->kappa - a->alpha) / a->K_N)) +
        (1. - a->K_ND) * exp(-(a->kappa - a->alpha) / a->K_N) *
        (1. - exp(-(a->lambda - a->kappa) / a->K_ND)) +
        a->mN * (a->lambda - a->kappa + (1. - a->K_ND) * (
            1. - exp(-(a->lambda - a->kappa) / a->K_ND)))));
    */

    a->Ixn_ = xy_pat[0] + n_pat[0] + d_pat[0];
    a->Ixy = xy_pat[1] + n_pat[1] + d_pat[1];
    a->Iyn = xy_pat[2] + n_pat[2] + d_pat[2];
    a->Iyd = xy_pat[3] + n_pat[3] + d_pat[3];
    a->Iynd = xy_pat[4] + n_pat[4] + d_pat[4];
    a->Iyn_ = xy_pat[5] + n_pat[5] + d_pat[5];
    a->Iyd_ = xy_pat[6] + n_pat[6] + d_pat[6];

    /**
     * Tree destructors
     */
    PopNode_free(X);
    PopNode_free(Y0);
    PopNode_free(Y1);
    PopNode_free(Y2);
    PopNode_free(N);
    PopNode_free(D);
    PopNode_free(ND);
    PopNode_free(XY);
    PopNode_free(AEND);
    gsl_rng_free(rng);
    return 0;
}

int main(void) {

    long double      twoN0 = 2.0*2308.0; /* haploid size of ancestral human pop */
    long double      gen = 29.0;         /* generation time */
    long double      s = twoN0 * gen;    /* years per time unit */
#if 0
    /* for production */
    int         nthreads = 12;   /* number of threads to launch */
    int         nTasks = 50;     /* total number of tasks */
    unsigned long   nreps = 10000000;
#else
    /* for debugging */
    int         nthreads = 4;   /* number of threads to launch */
    int         nTasks = 10;     /* total number of tasks */
    unsigned long   nreps = 10000;
#endif

    TaskArg targ = {
        .mN = 0.05,          /* primary (Neanderthal) admixture */
        .mD = 0.02,         /* secondary (Denisovan) admixture */

        /* Time backwards from the present, units of twoN0 gen */
        .alpha = 25e3/s,     /* Denisovan admixture */
        .delta = 55e3/s,     /* Neanderthal admixture */
        .zeta = 110e3/s,     /* X and Y split */
        .kappa = 427e3/s,    /* N and D split */
        .lambda = 658e3/s,   /* archaics and moderns split */
	    .lambda_max = 2 * 658e3/s, /* max value of lambda for
                                      perturbations */

        .mutRate = 1.48e-8 * twoN0, // mutation rate

        /* population sizes relative to N0 */
        .K_D = 1.0,
        .K_N = 1.0,
        .K_ND = 1.0,
        .K_XY = 1.0,

        .rng_seed = 0,
        .nreps = nreps,

        /* Returned values */
	    .Ixn=0.0,
        .oxn=0.0,
	    .Ixnd=0.0,
        .oxnd=0.0,
	    .Ixn_=0.0,
        .oxn_=0.0,
	    .Iyn=0.0,
        .oyn=0.0,
	    .Iynd=0.0,
        .oynd=0.0,
	    .Iyn_=0.0,
        .oyn_=0.0,
	    .Ixy=0.0,
        .oxy=0.0
    };

#if 1
    /* Populations differ in size */
    targ.K_D = 0.1;
    targ.K_N = 0.25;
    targ.K_ND = 0.25;
    targ.K_XY = 6.5;
#elif 0
    /* Force coalescent events into ancestral population */
    targ.K_D = targ.K_IMND = targ.K_IMN = strtod("Inf", 0);
#endif

    TaskArg *taskarg[nTasks];
    int         j, k;
    time_t      currtime = time(NULL);

    if(nthreads == 0)
        nthreads = getNumCores();

    if(nthreads > nTasks)
        nthreads = nTasks;

    printf("simultaneous_likelihood\n");
    printf("%-20s: %Lf\n", "2N0", twoN0);
    printf("%-20s: %Lf\n", "generation", gen);
    printf("%-20s: %d\n", "nthreads", nthreads);
    printf("%-20s: %d\n", "nTasks", nTasks);
    TaskArg_printParams(&targ, stdout);

    for(j = 0; j < nTasks; ++j)
        taskarg[j] = TaskArg_new(&targ, (unsigned) currtime+j);

    JobQueue   *jq = JobQueue_new(nthreads);
    if(jq == NULL)
        eprintf("ERR@%s:%d: Bad return from JobQueue_new",
                __FILE__, __LINE__);
    for(j = 0; j < nTasks; ++j)
        JobQueue_addJob(jq, taskfun, taskarg[j]);
    JobQueue_waitOnJobs(jq);
    fflush(stdout);
    fprintf(stderr, "Back from threads\n");

    printf("%8s %8s %8s %8s %8s %8s %8s "
           "%8s %8s %8s %8s %8s %8s %8s "
           "%8s %8s %8s %8s %8s %8s %8s\n",
           "Ixn", "oxn",
           "Ixn_", "oxn_",
           "Ixy", "oxy",
           "Iyn", "oyn",
           "Iyd", "oyd",
           "Iynd", "oynd",
           "Ixnd", "oxnd",
           "Iyn_", "oyn_",
           "alpha", "delta", "zeta", "kappa", "lambda");
    for(j=0; j<nTasks; ++j) {
        printf("%8.4Lf %8.4Lf %8.4Lf %8.4Lf %8.4Lf "
               "%8.4Lf %8.4Lf %8.4Lf %8.4Lf %8.4Lf "
               "%8.4Lf %8.4Lf %8.4Lf %8.4Lf %8.4Lf "
               "%8.4Lf %8Lf %8Lf %8Lf %8Lf %8Lf",
			taskarg[j]->Ixn, taskarg[j]->oxn,
			taskarg[j]->Ixn_, taskarg[j]->oxn_,
			taskarg[j]->Ixy, taskarg[j]->oxy,
			taskarg[j]->Iyn, taskarg[j]->oyn,
			taskarg[j]->Iyd, taskarg[j]->oyd,
			taskarg[j]->Iynd, taskarg[j]->oynd,
			taskarg[j]->Ixnd, taskarg[j]->oxnd,
			taskarg[j]->Iyn_, taskarg[j]->oyn_,
			taskarg[j]->alpha,
			taskarg[j]->delta,
			taskarg[j]->zeta,
			taskarg[j]->kappa,
			taskarg[j]->lambda);
        putchar('\n');
    }

    printf("Simulated counts of site patterns\n");
    printf(" %8s %8s %8s %8s %8s %8s\n",
           "Ixnd", "Ixn", "Ixy", "Iyn", "Iyd", "Iynd"); 
    for(j=0; j<nTasks; ++j) {
        for(k=0; k<N_SITE_PAT; ++k)
            printf(" %8ld", taskarg[j]->sitePatCount[k]);
        putchar('\n');
    }

    double x[nTasks], y[nTasks];

    /* plot yn */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn;
        y[j] = taskarg[j]->oyn;
    }
    pictex(x, y, nTasks, "Iyn", "oyn", "$Iyn$", "fig_yn.tex");

    /**
     * Alpha plots
     */

    /* plot (Iyn - Oyn, delta - alpha) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->delta - taskarg[j]->alpha;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\delta - \\alpha$",
               "$\\text{Diff}(\\delta)$", "fig_yn-da.tex");

    /* plot (Iyn - Oyn, zeta - alpha) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->zeta - taskarg[j]->alpha;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\zeta - \\alpha$",
               "$\\text{Diff}(\\zeta)$", "fig_yn-za.tex");

    /* plot (Iyn - Oyn, kappa - alpha) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->kappa - taskarg[j]->alpha;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\kappa - \\alpha$",
               "$\\text{Diff}(\\kappa)$", "fig_yn-ka.tex");

    /* plot (Iyn - Oyn, lambda - alpha) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->lambda - taskarg[j]->alpha;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\lambda - \\alpha$",
               "$\\text{Diff}(\\lambda)$", "fig_yn-la.tex");

    /**
     * Delta plots
     */
    /* plot (Iyn - Oyn, zeta - delta) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->zeta - taskarg[j]->delta;
    }

    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\zeta - \\delta$",
               "$\\text{Diff}(\\zeta)$", "fig_yn-zd.tex");

    /* plot (Iyn - Oyn, kappa - delta) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->kappa - taskarg[j]->delta;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\kappa - \\delta$",
               "$\\text{Diff}(\\kappa)$", "fig_yn-kd.tex");

    /* plot (Iyn - Oyn, lambda - delta) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->lambda - taskarg[j]->delta;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\lambda - \\delta$",
               "$\\text{Diff}(\\lambda)$", "fig_yn-ld.tex");

    /**
     * Zeta plots
     */
    /* plot (Iyn - Oyn, kappa - zeta) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->kappa - taskarg[j]->zeta;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\kappa - \\zeta$",
               "$\\text{Diff}(\\kappa)$", "fig_yn-kz.tex");

    /* plot (Iyn - Oyn, lambda - zeta) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->lambda - taskarg[j]->zeta;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\lambda - \\zeta$",
               "$\\text{Diff}(\\lambda)$", "fig_yn-lz.tex");

    /**
     * Kappa plots
     */
    /* plot (Iyn - Oyn, lambda - kappa) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->lambda - taskarg[j]->kappa;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$\\lambda - \\kappa$",
               "$\\text{Diff}(\\lambda)$", "fig_yn-lk.tex");

    /**
     * Population plots
     */
    /* plot (Iyn - Oyn, K_N) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->K_N;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$K_N$", "$\\text{Diff}(\\lambda)$",
               "fig_yn-Kn.tex");

    /* plot (Iyn - Oyn, K_D) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->K_D;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$K_D$", "$\\text{Diff}(\\lambda)$",
               "fig_yn-Kd.tex");

    /* plot (Iyn - Oyn, K_ND) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->K_ND;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$K_{ND}$", "$\\text{Diff}(\\lambda)$",
               "fig_yn-Knd.tex");

    /* plot (Iyn - Oyn, K_XY) */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iyn - taskarg[j]->oyn;
        y[j] = taskarg[j]->K_XY;
    }
    pictex_par(x, y, nTasks, "Iyn-oyn", "$K_{XY}$", "$\\text{Diff}(\\lambda)$",
               "fig_yn-Kae.tex");


    /* Thread destructor */
    for(j=0; j<nTasks; ++j)
        TaskArg_free(taskarg[j]);
    JobQueue_free(jq);

    return 0;
}
