/**
 * @file boot.h
 * @author Alan R. Rogers
 * @brief Header for boot.c.
 * @copyright Copyright (c) 2016 Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LEGO_BOOT_H
#  define LEGO_BOOT_H

#  include "typedefs.h"
#  include <gsl/gsl_rng.h>

BootChr    *BootChr_new(long nsnp, long nrep, long blockLength,
                           gsl_rng * rng);

#  ifndef NDEBUG
void        BootChr_sanityCheck(const BootChr *boot, const char *file, int line);
#  endif
long        BootChr_multiplicity(const BootChr * boot, long ndx, long rep);
void        BootChr_add(BootChr * self, long snp, int pat, double z);
long        BootChr_nrep(const BootChr * boot);
int         BootChr_nbin(const BootChr * boot);
long        BootChr_nblock(const BootChr * boot);
long        BootChr_nsnp(const BootChr * boot);
void        BootChr_free(BootChr * boot);
void        BootChr_get_rep(BootChr * self, int npat, double count[npat],
                            int rep);
void        BootChr_print(const BootChr * boot, FILE * ofp);

double      interpolate(double p, double *v, long len);
void        confidenceBounds(double *lowBnd, double *highBnd,
                             double confidence, double *v, long len);

#  ifndef NDEBUG
unsigned    BootChr_multiplicity_slow(BootChr * boot, long snp, long rep);
#  endif

#if 0
int         BootChr_equals(const BootChr * x, const BootChr * y);
BootChr    *BootChr_dup(const BootChr * old);
void        BootChr_plus_equals(BootChr * x, const BootChr * y);
long unsigned BootChr_rawCounts(const BootChr * boot, int rep, int bin,
                                double *numerator, double *denominator,
                                double *sumRsq, double *sep_cm);
long        BootChr_purge(BootChr * boot);

BootConf   *BootConf_new(BootChr * boot, double confidence);
void        BootConf_printHdr(const BootConf * bc, FILE * ofp);
double      BootConf_lowBound(const BootConf * bc, long bin);
double      BootConf_highBound(const BootConf * bc, long bin);
double      BootConf_loSpecBound(const BootConf * bc, long i);
double      BootConf_hiSpecBound(const BootConf * bc, long k);
void        BootConf_print(const BootConf * bc, FILE * ofp);
void        BootConf_free(BootConf * bc);

#endif
#endif
