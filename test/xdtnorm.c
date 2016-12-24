//  Example for using rtnorm
//
// This code was translated by Alan R. Rogers into C, based on the C++
// original written by G. Dolle and V. Mazet, which is available at
// http://miv.u-strasbg.fr/mazet/rtnorm/rtnormCpp.zip
//
//  Copyright (C) 2012 Guillaume Dollé, Vincent Mazet (LSIIT,
//  CNRS/Université de Strasbourg)
//  Licence: GNU General Public License Version 2
//  see http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
//
//  Depends: LibGSL
//  OS: Unix based system

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "dtnorm.h"

typedef struct Pair Pair;
struct Pair {
    double first, second;
};

void test(int nreps, double mu, double sigma, double a, double b,
          gsl_rng *rng, int verbose);
Pair meanAndLeftHalf(double u, double s, double a, double b);

/**
   Return: .first = mean; .second = fraction of probability mass in
   left half of interval (a,b).

   alpha = (a-u)/sigma
   beta  = (b-u)/sigma
   f(x) = density
   F(x) = cumulative density
   Z = F(beta) - F(alpha)

   E[x] = mu + sigma*(f(alpha) - f(beta))/(F(beta) - F(alpha))

   gsl_ran_ugaussian_pdf(x)  // mean=zero, sigma=1
   gsl_cdf_ugaussian_P(x)
 */
Pair meanAndLeftHalf(double u, double s, double a, double b) {
    double alpha = (a-u)/s;
    double beta = (b-u)/s;
    double fa = gsl_ran_ugaussian_pdf(alpha);
    double fb = gsl_ran_ugaussian_pdf(beta);
    double Fa = gsl_cdf_ugaussian_P(alpha);
    double Fb = gsl_cdf_ugaussian_P(beta);
    double Fmid = gsl_cdf_ugaussian_P(0.5*(alpha+beta));
    double Z = Fb - Fa;
    double mean = (fa - fb)/Z;

    Pair p = {
        .first = u + s*mean,
        .second = (Fmid-Fa)/Z
    };
    return p;
}

void test(int nreps, double mu, double sigma, double a, double b,
          gsl_rng *rng, int verbose) {
    if(verbose)
        printf("test(%lf, %lf, %lf, %lf)\n", mu, sigma, a, b);
    int i, nlefthalf=0;
    double min = HUGE_VAL;
    double max = -HUGE_VAL;
    double x, m=0, lefthalf=0;
    double mid = 0.5*(a+b);
    for(i=0; i < nreps; ++i) {
        x = dtnorm(mu, sigma, a, b, rng);
        m += x;
        if(x < mid)
            ++nlefthalf;
        min = fmin(min, x);
        max = fmax(max, x);
    }
    assert(min >= a);
    assert(max <= b);
    m /= nreps;
    lefthalf = nlefthalf / ((double) nreps);
    Pair p = meanAndLeftHalf(mu, sigma, a, b);
    if(verbose) {
        printf("   mean : obs=%lf expected=%lf diff=%lf\n",
               m, p.first, fabs(m-p.first));
        printf("   lhalf: obs=%lf expected=%lf diff=%lf\n",
               lefthalf, p.second, fabs(lefthalf-p.second));
    }
    assert(fabs(p.first - m) < 0.005*fmax(1.0, p.first));
    assert(fabs(p.second - lefthalf) < 0.005*fmax(1.0, p.second));
    if(verbose)
        printf("   Passed\n");
}

int main(int argc, char **argv) {
    int         verbose = 0;
    double      a;              // Left bound
    double      b;              // Right bound
    double      mu = 0.1;       // Mean
    double      sigma = 0.75;   // Standard deviation
    int         nreps = 1e6;    // Number of random variates to generate

    if(argc==2 && 0==strcmp("-v", argv[1]))
        verbose = 1;
    else if(argc != 1) {
        fprintf(stderr,"usage: xdtnorm [-v]\n");
        exit(1);
    }

    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned long) time(NULL));

    // Test exponential proposal
    a = 3.1;
    b = a + 0.75;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);

    // Test gaussian proposal
    a = -1.0;
    b = a + 1.5;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);

    // Test Robert's method
    a = -1.0;
    b = -0.5;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);
    a = -0.5;
    b = 0.2;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);
    a = 3.1;
    b = a + 0.3;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);

    // Test Chopin's method
    a = 0.1;
    b = a + 0.05;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);
    a = 0.1;
    b = a + 2.0;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);
    a = 2.8;
    b = a + 0.05;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);
    a = 2.8;
    b = a + 2.0;
    test(nreps, mu, sigma, mu+a*sigma, mu+b*sigma, rng, verbose);

    printf("dtnorm OK\n");

    gsl_rng_free(rng);

    return 0;
}
