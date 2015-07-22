/**
 * @file q.c
 * @brief Calculate site pattern frequencies under model for Q.
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "misc.h"
#include "gptree.h"
#include <stdio.h>

void segment(const char *name, double t, double n, int samples);
void mix(const char *ego, const char *ancestor, double m,
         const char *introgressor);
void derive(const char *ego, const char *ancestor);


void segment(const char *name, double t, double n, int samples) {
    printf("segment %3s t=%9.6lf N=%9.6lf", name, t, n);
    if(samples)
        printf(" samples=%d", samples);
    putchar('\n');
}

void mix(const char *ego, const char *ancestor, double m,
         const char *introgressor) {
    printf("mix %s from %9.6lf %s + %9.6lf %s\n", ego, 1-m, ancestor, m,
           introgressor);
}

void derive(const char *ego, const char *ancestor) {
    printf("derive %s from %s\n", ego, ancestor);
}

int main(void) {

    double      twoN0 = 2e4;     // haploid size of ancestral human pop
    double      gen = 29.0;      // generation time
    double      s = twoN0 * gen; // years per time unit


    double mN = 0.02;            // N->EV gene flow
    double mD = 0.01;            // D->V gene flow

    // Time backwards from the present, units of twoN0 gen
    double alpha = 25e3/s;       // Denisovan admixture
    double delta = 55e3/s;       // Neanderthal admixture
    double epsilon = 65e3/s;     // age of older Neanderthal fossil
    double zeta = 110e3/s;       // Y and EV split
    double kappa = 427e3/s;      // N and D split
    double lambda = 658e3/s;     // archaics and moderns split

    // population sizes relative to N0
    double K_X, K_Y, K_D, K_XY, K_ND, K_N;
    K_X = K_Y = K_D = K_XY = K_ND = K_N = 1.0;

#if 1
    // Populations differ in size
    K_X = 2.0;
    K_XY = 2.0;
    K_ND = 0.2;
    K_N = 0.1;
    K_D = 0.1;
#elif 0
    // Force coalescent events into ancestral population
    K_X = K_Y = K_D = K_XY = K_ND = K_N = strtod("Inf", 0);
#endif

    // survival probabilities (s) and their complements (f)
    double sND = survival(lambda - kappa, K_ND);
    double fND = 1.0 - sND;
    double sN  = survival(kappa - delta, K_N);
    double fN  = 1.0 - sN;
    double sXY = survival(lambda - zeta, K_XY);
    double fXY = 1.0 - sXY;

    // Expected fraction of each site pattern
    double Inx = (mN*(1.0-mD)*sN*sND + mD*sND + (1.0 - mN)*(1.0 - mD)*sXY)/3.0;
    double Iny = Inx + mN*(1.0 - mD)*(lambda - delta + (1.0 - K_N)*fN
                                      + (1.0 - K_ND)*sN*fND)
        + mD*(lambda - kappa + (1.0 - K_ND)*fND);
    double Ixy = Inx + (1.0 - mN)*(1.0 - mD)*(lambda - zeta + (1.0 - K_XY)*fXY);

    double sum = Inx+Iny+Ixy;
    Inx /= sum;
    Iny /= sum;
    Ixy /= sum;

    printf("# Q\n");
    printf("# 2N0         : %0.0lf\n", twoN0);
    printf("# t/2N0       : delta=%lf alpha=%lf epsilon=%lf\n", delta, alpha, epsilon);
    printf("#             : kappa=%lf lambda=%lf zeta=%lf\n", kappa, lambda, zeta);
    printf("# Rel pop size: X=%lf Y=%lf N=%lf D=%lf\n", K_X, K_Y, K_N, K_D);
    printf("#             : XY=%lf ND=%lf URHUMAN=%lf\n", K_XY, K_ND, 1.0);
    printf("# admixture->Y: mN=%lf mD=%lf\n", mN, mD);
    printf("# Inx         : %lf\n", Inx);
    printf("# Iny         : %lf\n", Iny);
    printf("# Ixy         : %lf\n", Ixy);

    segment("x", 0, K_X, 1);
    segment("y", 0, K_Y, 1);
    segment("w", alpha, K_Y, 0);
    segment("z", delta, K_Y, 0);
    segment("n", delta, K_N, 1);
    segment("d", alpha, K_D, 1);
    segment("xy", zeta, K_XY,0);
    segment("nd", kappa, K_ND,0);
    segment("a", lambda, 1.0, 0);
    mix("y", "w", mD, "d");
    mix("w", "z", mN, "n");
    derive("x", "xy");
    derive("z", "xy");
    derive("xy", "a");
    derive("n", "nd");
    derive("d", "nd");
    derive("nd", "a");
    
    return 0;
}
