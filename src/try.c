#include <assert.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "rtnorm.h"

// Code for Chopin's algorithm is at
// https://sites.google.com/site/nicolaschopinstatistician/software 

/*
  Prefer Robert if b+a <= 0.0 or if (b-a < 1 and b+a < 2)
  or if b-a < 0.25

b = (bpa + bma)/2
bpa - (bpa + bma)/2 = k = (bpa - bma)/2
2*k = bpa - bma
bma = -2*k + bpa

*/

double timefun(int nreps, double a, double b,
               double (*fun)(double m, double sd, double a, double b, gsl_rng *rng),
               gsl_rng *rng);
double dtnorm_simple(double mean, double sd, double a, double b, gsl_rng *rng);
double dtnorm_robert(double mean, double sd, double a, double b, gsl_rng *rng);
double dtnorm(double mean, double sd, double a, double b, gsl_rng *rng);

double dtnorm(double mean, double sd, double a, double b, gsl_rng *rng) {
    if(a >= b) {
        fprintf(stderr, "%s:%d: b <= a\n", __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    if( (b-a)/sd < 1.25 )
        return dtnorm_robert(mean, sd, a, b, rng);
    return dtnorm_simple(mean, sd, a, b, rng);
}

double dtnorm_simple(double mean, double sd, double a, double b, gsl_rng *rng) {

    double x;
    double da = a - mean;
    double db = b - mean;

    do{
        x = gsl_ran_gaussian_ziggurat(rng, sd);
    }while(x < da || x > db);

    return mean + x;
}

/// Sample from a doubly-truncated normal distribution, with support
/// on [a,b].  Parameters mean and sd refer to the underlying
/// untruncated normal distribution.  Uses algorithm of Robert,
/// CP. 1995. Simulation of truncated normal variables. Statistics and
/// Computing, 5(2):121-125.
double dtnorm_robert(double mean, double sd, double a, double b, gsl_rng *rng) {
    if(b <= a) {
        fprintf(stderr, "%s:%d: b <= a\n", __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    double za = (a - mean)/sd;
    double zb = (b - mean)/sd;
    double u, z, rho;

    do{
        z = gsl_ran_flat(rng, za, zb);
        if( za*zb < 0.0 )               // of opposite sign
            rho = exp(-0.5*z*z);
        else if( zb < 0.0 )             // both negative
            rho = exp( 0.5*(zb*zb - z*z) );
        else{                           // both positive
            assert(za > 0.0);
            rho = exp( 0.5*(za*za - z*z) );
        }

        u = gsl_rng_uniform(rng);
    }while( u > rho );

    return mean + z*sd;
}

double timefun(int nreps, double a, double b,
               double (*fun)(double m, double sd, double a, double b, gsl_rng *rng),
               gsl_rng *rng) {
    clock_t start, finish;
    int i;

    start = clock();
    for(i=0; i < nreps; ++i)
        (void) (*fun)(0.0, 1.0, a, b, rng);
    finish = clock();
    return ((double) (finish-start)) / CLOCKS_PER_SEC;
}

int main(int argc, char **argv) {

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned long) time(NULL));
    
    int trial, ntrials=1000, nreps=10000000;
    double a, b, bma, bpa;
    double trobert, tchopin;
    //    double amin=1.5, amax = 3.0, diffmax = 3.0;
    double bpa0=-5.0, bpa1 = 10.0;
    double bma0 = 0.01, bma1 = 5.0;

    printf("%7s %7s %7s %7s\n", "a", "b", "tr", "tc");
    for(trial = 0; trial < ntrials; ++trial) {
        bpa = gsl_ran_flat(rng, bpa0, bpa1);       // b + a
        bma = gsl_ran_flat(rng, bma0, bma1);       // b - a
        b = (bpa + bma)/2.0;
        a = bpa - b;
        
        trobert = timefun(nreps, a, b, rtnorm_robert, rng);
        tchopin = timefun(nreps, a, b, rtnorm, rng);

        printf("%7.4lf %7.4lf %7.4lf %7.4lf\n", a, b, trobert, tchopin);
        fflush(stdout);
    }

    gsl_rng_free(rng);
	return 0;
}
