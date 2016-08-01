/**
 * @file boot.h
 * @author Alan R. Rogers
 * @brief Header for boot.c.
 * @copyright Copyright (c) 2016 Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LEGO_BOOT_H
#define LEGO_BOOT_H

#include "typedefs.h"
#include <gsl/gsl_rng.h>
void        confidenceBounds(double *lowBnd, double *highBnd,
                             double confidence, double *v, long len);
double      interpolate(double p, double *v, long len);
Boot       *Boot_new(long nSNPs, long nReps, unsigned twoNsmp,
                     int folded, long blockLength,
                     double windowcm, int nBins, gsl_rng * rng);

// Add contribution z to i'th site pattern
void        Boot_add(Boot * boot, int i, double z);

void        Boot_free(Boot * boot);
int         Boot_equals(const Boot * x, const Boot * y);
Boot       *Boot_dup(const Boot * old);
long        Boot_nReps(const Boot * boot);
int         Boot_nBins(const Boot * boot);
long        Boot_nBlocks(const Boot * boot);
long        Boot_nSNPs(const Boot * boot);
void        Boot_plus_equals(Boot * x, const Boot * y);
#ifndef NDEBUG
void        Boot_sanityCheck(const Boot * boot, const char *file, int line);
#endif
long        Boot_multiplicity(const Boot * boot, long ndx, long rep);
void        Boot_get_rep(Boot * boot, DblArray *sigdsq, DblArray *rsq,
                         DblArray *cm, ULIntArray *nobs,
                         ULIntArray *spectrum, int rep);
long unsigned Boot_rawCounts(const Boot * boot, int rep, int bin,
                             double *numerator, double *denominator,
                             double *sumRsq, double *sep_cm);
long        Boot_purge(Boot * boot);
void        Boot_print(const Boot * boot, FILE * ofp);

#ifndef NDEBUG
unsigned Boot_multiplicity_slow(Boot * boot, long snp, long rep);
#endif

BootConf   *BootConf_new(Boot * boot, double confidence);
void        BootConf_printHdr(const BootConf * bc, FILE * ofp);
double      BootConf_lowBound(const BootConf * bc, long bin);
double      BootConf_highBound(const BootConf * bc, long bin);
double      BootConf_loSpecBound(const BootConf * bc, long i);
double      BootConf_hiSpecBound(const BootConf * bc, long k);
void        BootConf_print(const BootConf * bc, FILE * ofp);
void        BootConf_free(BootConf * bc);

#endif
