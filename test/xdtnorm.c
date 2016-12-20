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
#include <time.h>

#include "dtnorm.h"

int main(void) {
    double      a;              // Left bound
    double      b;              // Right bound
    double      mu = 2;         // Mean
    double      sigma = 3;      // Standard deviation
    double      x;              // Output of rtnorm
    int         nreps = 1e6;    // Number of random variates to generate

    //--- GSL random init ---
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned long) time(NULL));

    //--- generate and display the random numbers ---
    printf("underlying distribution: Normal(%lf, %lf)\n", mu, sigma);
    for(int k = 0; k < nreps; k++) {
        a = gsl_ran_flat(rng, mu - 5.0*sigma, mu + 5.0*sigma);
        b = gsl_ran_flat(rng, a + 0.001*sigma, a + 5.0*sigma);
        x = dtnorm(mu, sigma, a, b, rng);
        printf("%lf\n", x);
    }

    gsl_rng_free(rng);

    return 0;
}
