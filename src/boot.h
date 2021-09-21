/**
 * @file boot.h
 * @author Alan R. Rogers
 * @brief Header for boot.c.
 * @copyright Copyright (c) 2016,2018 Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef LEGO_BOOT_H
#  define LEGO_BOOT_H

#  include "typedefs.h"
#  include <gsl/gsl_rng.h>

#define DEF_BOOT_BLOCK_SIZE 500

Boot       *Boot_new(int nchr, long nsnp[nchr], long nrep, int npat,
                     long blocksize, gsl_rng * rng);
void        Boot_free(Boot * self);
void        Boot_add(Boot *self, int chr, long snpndx, int pat, double z);
void        Boot_aggregate(Boot * self, int rep, int npat,
                           double count[npat]);
void        Boot_print(const Boot * self, FILE * ofp);
#ifndef NDEBUG
void        Boot_sanityCheck(const Boot * self, const char *file, int line);
#endif

void        confidenceBounds(double *lowBnd, double *highBnd, double confidence,
							 long len, double v[len]);

#  ifndef NDEBUG
unsigned    Boot_multiplicity_slow(Boot * self, long snp, long rep);
#  endif
#endif
