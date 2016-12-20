//  Pseudorandom numbers from a truncated Gaussian distribution.
//
// This code was translated by Alan R. Rogers into C, based on the C++
// original written by G. Dolle and V Mazet, which is available at 
// http://miv.u-strasbg.fr/mazet/rtnorm/rtnormCpp.zip
//
//  This implements an extension of Chopin's algorithm detailed in
//  N. Chopin, "Fast simulation of truncated Gaussian distributions",
//  Stat Comput (2011) 21:275-288
//  
//  Copyright (C) 2012 Guillaume Dollé, Vincent Mazet
//  (LSIIT, CNRS/Université de Strasbourg)
//  Version 2012-07-04, Contact: vincent.mazet@unistra.fr
//  
//  06/07/2012:
//  - first launch of rtnorm.cpp
//  
//  Licence: GNU General Public License Version 2
//  This program is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2 of the License, or (at your
//  option) any later version. This program is distributed in the hope that
//  it will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details. You should have received a
//  copy of the GNU General Public License along with this program; if not,
//  see http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
//  
//  Depends: LibGSL
//  OS: Unix based system

#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include <assert.h>

#include "dtnorm.h"
#include "dtnorm_data.h"

int         N = 4001;           // Index of the right tail

// Pseudorandom numbers from a truncated Gaussian distribution
// The Gaussian has parameters mu (default 0) and sigma (default 1)
// and is truncated on the interval [a,b].
// Returns the random variable x and its probability p(x).
double dtnorm(const double mu, const double sigma, double a, double b, gsl_rng *rng) {
    // Design variables
    //double      xmin = -2.00443204036;  // Left bound
    double      xmax = 3.48672170399;   // Right bound
    int         kmin = 5;       // if kb-ka < kmin then use a rejection algorithm
    double      INVH = 1631.73284006;   // = 1/h, h being the minimal interval range
    int         I0 = 3271;      // = - floor(x(0)/h)
    double      ALPHA = 1.837877066409345;  // = log(2*pi)
    int         xsize = sizeof(x) / sizeof(double); // Length of table x
    int         stop = false;

    double      r, z, e, ylk, simy, lbound, u, d, sim;
    int         i, ka, kb, k;

    // Scaling
    if(mu != 0 || sigma != 1) {
        a = (a - mu) / sigma;
        b = (b - mu) / sigma;
    }

    double bma = b - a;

    // Check if a < b
    if(a >= b) {
        fprintf(stderr, "%s:%d: *** B must be greater than A ! ***\n",
                __FILE__, __LINE__);
        exit(1);
    }

    // Check if |a| < |b|
    else if(fabs(a) > fabs(b))
        r = -dtnorm(0.0, 1.0, -b, -a, rng);

    // Truncated exponential proposal   
    else if(a > 3.0 && bma > 0.5)
        r = dtexp(rng, a, b);

    // Gaussian proposal 
    else if(a < 0.0 && b > 0.0 && bma > 1.0) {
        while(!stop) {
            r = gsl_ran_gaussian_ziggurat(rng, 1.0);
            stop = (r >= a) && (r <= b);
        }
    }

    // Robert's algorithm
    else if(b < 0.0
            || (a < 0.0 && (b < 0.0 || bma < 1.0))
            || (a > 3.0 && bma < 0.5)) {
        double rho;
        do{
            r = gsl_ran_flat(rng, a, b);
            if( a*b < 0.0 )               // of opposite sign
                rho = exp(-0.5*r*r);
            else if( b < 0.0 )             // both negative
                rho = exp( 0.5*(b*b - r*r) );
            else{                           // both positive
                assert(a > 0.0);
                rho = exp( 0.5*(a*a - r*r) );
            }

            u = gsl_rng_uniform(rng);
        }while( u > rho );
    }

    // Chopin's algorithm
    else {
        if(!(0.0 <= a && a <= 3.0)) {
            fprintf(stderr,"%s:%d: (a,b) = (%lf, %lf)\n",
                    __FILE__, __LINE__, a, b);
            exit(1);
        }
        assert(0.0 <= a);
        assert(a <= 3.0);
        // Compute ka
        i = I0 + floor(a * INVH);
        ka = ncell[i];

        // Compute kb
        (b >= xmax) ?
            kb = N : (i = I0 + floor(b * INVH), kb = ncell[i]
            );

        // If |b-a| is small, use rejection algorithm with a truncated
        // exponential proposal 
        if(abs(kb - ka) < kmin) {
            r = dtexp(rng, a, b);
            stop = true;
        }

        while(!stop) {
            // Sample integer between ka and kb
            k = floor(gsl_rng_uniform(rng) * (kb - ka + 1)) + ka;

            if(k == N) {
                // Right tail
                lbound = x[xsize - 1];
                z = -log(gsl_rng_uniform(rng));
                e = -log(gsl_rng_uniform(rng));
                z = z / lbound;

                if((z*z <= 2 * e) && (z < b - lbound)) {
                    // Accept this proposition, otherwise reject
                    r = lbound + z;
                    stop = true;
                }
            }

            else if((k <= ka + 1) || (k >= kb - 1 && b < xmax)) {

                // Two leftmost and rightmost regions
                sim = x[k] + (x[k + 1] - x[k]) * gsl_rng_uniform(rng);

                if((sim >= a) && (sim <= b)) {
                    // Accept this proposition, otherwise reject
                    simy = yu[k] * gsl_rng_uniform(rng);
                    if((simy < yl(k))
                       || (sim * sim + 2 * log(simy) + ALPHA) < 0) {
                        r = sim;
                        stop = true;
                    }
                }
            }

            else                // All the other boxes
            {
                u = gsl_rng_uniform(rng);
                simy = yu[k] * u;
                d = x[k + 1] - x[k];
                ylk = yl(k);
                if(simy < ylk)  // That's what happens most of the time 
                {
                    r = x[k] + u * d * yu[k] / ylk;
                    stop = true;
                } else {
                    sim = x[k] + d * gsl_rng_uniform(rng);

                    // Otherwise, check you're below the pdf curve
                    if((sim * sim + 2 * log(simy) + ALPHA) < 0) {
                        r = sim;
                        stop = true;
                    }
                }

            }
        }
    }

    // Scaling
    if(mu != 0 || sigma != 1)
        r = r * sigma + mu;

    return r;
}

// Compute y_l from y_k
double yl(int k) {
    double      yl0 = 0.053513975472;   // y_l of the leftmost rectangle
    double      ylN = 0.000914116389555;    // y_l of the rightmost rectangle

    if(k == 0)
        return yl0;

    else if(k == N - 1)
        return ylN;

    else if(k <= 1953)
        return yu[k - 1];

    else
        return yu[k + 1];
}

// Rejection algorithm with a truncated exponential proposal
double dtexp(gsl_rng * rng, double a, double b) {
    double      twoasq = 2*a*a;
    double      expab = expm1(-a * (b - a));
    double      z, e;

    do{
        z = log(1 + gsl_rng_uniform(rng) * expab);
        e = -log(gsl_rng_uniform(rng));
    }while(twoasq*e <= z*z);
    return a - z/a;
}
