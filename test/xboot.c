/**
 * @file xboot.c
 * @author Alan R. Rogers
 * @brief Test boot.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "boot.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <float.h>
#include <time.h>

int main(int argc, char **argv) {
    long        nsnp[] = {8, 10, 15};
    int         nchr = (sizeof nsnp)/(sizeof nsnp[0]);
    long        nReps = 10;
    int         npat = 10;
    long        blockLength = 3;
    int         verbose = 0;

    long        i, j;
    time_t      currtime = time(NULL);
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr,"usage: xboot [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr,"usage: xboot [-v]\n");
        exit(EXIT_FAILURE);
    }

    gsl_rng_set(rng, (unsigned) currtime);

    double      v[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };

    assert(Dbl_near(interpolate(0.0, v, 6), 0.0));
    assert(Dbl_near(interpolate(1.0, v, 6), 5.0));
    assert(Dbl_near(interpolate(0.5, v, 6), 2.5));
    assert(Dbl_near(interpolate(0.5, v, 5), 2.0));
    unitTstResult("interpolate", "OK");

    bootchr = BootChr_new(0, nsnp[0], nReps, npat, blockLength, rng);
    if(verbose)
        BootChr_print(bootchr, stdout);

    assert(nsnp[0] == BootChr_nsnp(bootchr));
    assert(nReps == BootChr_nrep(bootchr));
    assert(npat == BootChr_npat(bootchr));
    assert(BootChr_nblock(bootchr) ==
           floor(0.5+(nsnp[0]/((double) blockLength))));

    long        isnp, irep;

    for(isnp = 0; isnp < nsnp[0]; ++isnp) {
        for(irep = 0; irep < nReps; ++irep) {
            long        m1 = BootChr_multiplicity_slow(bootchr, isnp, irep);
            long        m2 = BootChr_multiplicity(bootchr, isnp, irep);

            assert(m1 == m2);
        }
    }
    unitTstResult("BootChr_new", "OK");

    for(i=0; i < npat; ++i) {
        long snp = gsl_rng_uniform_int(rng, nsnp[0]);
        BootChr_add(bootchr, snp, i, 1.0);
    }

    for(i=0; i < nReps; ++i) {
        double count[npat];
        memset(count, 0, sizeof count);
        BootChr_aggregate(bootchr, i, npat, count);
        if(verbose) {
            printf("rep %2ld:", i);
            for(j=0; j < npat; ++j)
                printf(" %6.2lf", count[j]);
            putchar('\n');
        }
    }

    if(verbose && nReps <= 5 && BootChr_nblock(bootchr) <= 50)
        BootChr_print(bootchr, stdout);

    long        snp, rep, m, slow;

    for(i = 0; i < 100; ++i) {
        rep = gsl_rng_uniform_int(rng, nReps);
        snp = gsl_rng_uniform_int(rng, nsnp[0]);
        m = BootChr_multiplicity(bootchr, snp, rep);
        slow = BootChr_multiplicity_slow(bootchr, snp, rep);
        if(m != slow)
            eprintf("Boot_multiplicity FAILED@$s:%d:"
                    "rep=%ld snp=%ld m=%ld != slow=%ld\n",
                    __FILE__, __LINE__, rep, snp, m, slow);
    }

    BootChr_free(bootchr);

    unitTstResult("BootChr", "OK");

    Boot *boot = Boot_new(nchr, nsnp, nReps, npat, blockLength, rng);

    for(i=0; i < nchr; ++i) {
        for(j=0; j < npat; ++j) {
            snp = gsl_rng_uniform_int(rng, nsnp[i]);
            Boot_add(boot, i, snp, j, 1.0);
        }
    }

    for(i=0; i < nReps; ++i) {
        double count[npat];
        memset(count, 0, sizeof count);
        Boot_aggregate(boot, i, npat, count);
        if(verbose) {
            printf("rep %2ld:", i);
            for(j=0; j < npat; ++j)
                printf(" %6.2lf", count[j]);
            putchar('\n');
        }
    }
    Boot_free(boot);

    unitTstResult("Boot", "OK");

    gsl_rng_free(rng);

    return 0;
}
