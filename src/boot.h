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

BootChr    *BootChr_new(long nsnp, long nrep, int npat, long blocksize,
                        gsl_rng * rng);

#  ifndef NDEBUG
void        BootChr_sanityCheck(const BootChr * self, const char *file,
                                int line);
#  endif
long        BootChr_multiplicity(const BootChr * self, long snpndx, long rep);
void        BootChr_add(BootChr * self, long snpndx, int pat, double z);
long        BootChr_nrep(const BootChr * self);
long        BootChr_nblock(const BootChr * self);
long        BootChr_nsnp(const BootChr * self);
long        BootChr_npat(const BootChr * self);
void        BootChr_free(BootChr * self);
void        BootChr_print(const BootChr * self, FILE * ofp);
void        BootChr_aggregate(BootChr * self, int rep, int npat,
                              double count[npat]);

Boot       *Boot_new(int nchr, long nsnp[nchr], long nrep, int npat,
                     long blocksize, gsl_rng * rng);
void        Boot_free(Boot * self);
void        Boot_add(Boot *self, int chr, long snpndx, int pat, double z);
void        Boot_aggregate(Boot * self, int rep, int npat,
                           double count[npat]);
#ifndef NDEBUG
void        Boot_sanityCheck(const Boot * self, const char *file, int line);
#endif

double      interpolate(double p, double *v, long len);
void        confidenceBounds(double *lowBnd, double *highBnd, double confidence,
							 long len, double v[len]);
long        adjustBlockLength(long lengthWanted, int nsnp);

#  ifndef NDEBUG
unsigned    BootChr_multiplicity_slow(BootChr * self, long snp, long rep);
#  endif

#  if 0
int         BootChr_equals(const BootChr * x, const BootChr * y);
BootChr    *BootChr_dup(const BootChr * old);
void        BootChr_plus_equals(BootChr * x, const BootChr * y);
long unsigned BootChr_rawCounts(const BootChr * self, int rep, int bin,
                                double *numerator, double *denominator,
                                double *sumRsq, double *sep_cm);
long        BootChr_purge(BootChr * self);

BootConf   *BootConf_new(BootChr * bootchr, double confidence);
void        BootConf_printHdr(const BootConf * bc, FILE * ofp);
double      BootConf_lowBound(const BootConf * bc, long bin);
double      BootConf_highBound(const BootConf * bc, long bin);
void        BootConf_print(const BootConf * bc, FILE * ofp);
void        BootConf_free(BootConf * bc);

#  endif
#endif
