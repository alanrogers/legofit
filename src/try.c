#include <assert.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include "rtn.h"

// Code for Chopin's algorithm is at
// https://sites.google.com/site/nicolaschopinstatistician/software 

/*
d <- read.table("tnorm.txt", header=T)
ggplot(d, aes(a, b-a, color=best)) + geom_point()

Exp if 3 < a  and b-a > 0.5

Gauss if a < 0 < b and b-a > 0.75

Robert if (a < 0 and b < 0 or b-a < 0.75) or (a > 3 and b-a < 0.5)

Chopin if 0 < a < 3

*/
double timefun(int nreps, double a, double b,
               double (*fun)(double a, double b, gsl_rng *rng),
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
               double (*fun)(double a, double b, gsl_rng *rng),
               gsl_rng *rng) {
    clock_t start, finish;
    int i;

    start = clock();
    for(i=0; i < nreps; ++i)
        (void) (*fun)(a, b, rng);
    finish = clock();
    return ((double) (finish-start)) / CLOCKS_PER_SEC;
}

int main(int argc, char **argv) {

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned long) time(NULL));
    
    double      xmin = -2.00443204036;  // Left bound
    double      xmax = 3.48672170399;   // Right bound
    int trial, ntrials=1000, nreps=10000000;
    double a, b, bma, bpa;
    double trobert, tchopin, tgauss, texp;
    //    double amin=1.5, amax = 3.0, diffmax = 3.0;
    double bpa0=-10.0, bpa1 = 10.0;
    double bma0 = 0.01, bma1 = 7.0;

    printf("%7s %7s %s\n", "a", "b", "best");
    for(trial = 0; trial < ntrials; ++trial) {
#if 0
        bpa = gsl_ran_flat(rng, bpa0, bpa1);       // b + a
        bma = gsl_ran_flat(rng, bma0, bma1);       // b - a
        b = (bpa + bma)/2.0;
        a = bpa - b;
#else
        a = gsl_ran_flat(rng, -8.0, 5.0);
        b = gsl_ran_flat(rng, a + 0.01, 8.0);
#endif
        
        trobert = timefun(nreps, a, b, robert, rng);

        if(a >= xmin && a <= xmax)
            tchopin = timefun(nreps, a, b, chopin, rng);
        else
            tchopin = HUGE_VAL;

        if(a < 0.0 && b-a > 0.5 && b > 0.0) {
            tgauss = timefun(nreps, a, b, gaussproposal, rng);
        }else
            tgauss = HUGE_VAL;

        if(a > 3.0) {
            texp = timefun(nreps, a, b, rtexp, rng);
        }else
            texp = HUGE_VAL;

        double best = fmin(trobert,tchopin);
        best = fmin(best, tgauss);
        best = fmin(best, texp);

        char buff[20];

        buff[0] = '\0';
        if(trobert == best)
            strcat(buff,"r");
        if(tchopin == best)
            strcat(buff,"c");
        if(tgauss == best)
            strcat(buff,"g");
        if(texp == best)
            strcat(buff,"e");

        printf("%7.4lf %7.4lf %s\n", a, b, buff);
        fflush(stdout);
    }

    gsl_rng_free(rng);
	return 0;
}
