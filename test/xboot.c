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

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <float.h>
#include <time.h>

int main(int argc, char **argv) {
    long        nSNPs = 8;
    long        nReps = 10;
    int         npat = 10;
    long        blockLength = 3;
    int         verbose = 0;

    long        i;
    time_t      currtime = time(NULL);
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    BootChr       *bootchr;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xboot [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xboot [-v]\n");
    }

    gsl_rng_set(rng, (unsigned) currtime);

    double      v[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };

    assert(Dbl_near(interpolate(0.0, v, 6), 0.0));
    assert(Dbl_near(interpolate(1.0, v, 6), 5.0));
    assert(Dbl_near(interpolate(0.5, v, 6), 2.5));
    assert(Dbl_near(interpolate(0.5, v, 5), 2.0));
    unitTstResult("interpolate", "OK");

    bootchr = BootChr_new(nSNPs, nReps, npat, blockLength, rng);
    if(verbose)
        BootChr_print(bootchr, stdout);

    long        isnp, irep;

    for(isnp = 0; isnp < nSNPs; ++isnp) {
        for(irep = 0; irep < nReps; ++irep) {
            long        m1 = BootChr_multiplicity_slow(bootchr, isnp, irep);
            long        m2 = BootChr_multiplicity(bootchr, isnp, irep);

            assert(m1 == m2);
        }
    }
    unitTstResult("BootChr_new", "OK");

    if(verbose && nReps <= 5 && BootChr_nblock(bootchr) <= 50)
        BootChr_print(bootchr, stdout);

    long        snp, rep, m, slow;

    for(i = 0; i < 100; ++i) {
        rep = gsl_rng_uniform_int(rng, nReps);
        snp = gsl_rng_uniform_int(rng, nSNPs);
        m = BootChr_multiplicity(bootchr, snp, rep);
        slow = BootChr_multiplicity_slow(bootchr, snp, rep);
        if(m != slow)
            eprintf("Boot_multiplicity FAILED@$s:%d:"
                    "rep=%ld snp=%ld m=%ld != slow=%ld\n",
                    __FILE__, __LINE__, rep, snp, m, slow);
    }
    unitTstResult("BootChr_multiplicity", "OK");

#if 0
    double      s;
    long        ndx1, ndx2;
    int         polymorphic = 0;
    unsigned    gtypeSize = 10;
    SNP        *snp1 = SNP_new(gtypeSize, BootChr_nReps(boot));
    SNP        *snp2 = SNP_new(gtypeSize, BootChr_nReps(boot));
    unsigned char gtype[gtypeSize];
    int         nGtype;

    for(i = 0; i < 100; ++i) {
        ndx1 = gsl_rng_uniform_int(rng, nSNPs);
        ndx2 = gsl_rng_uniform_int(rng, nSNPs);
        nGtype = encodeHaploid(gtype, sizeof(gtype), "0101010101");
        assert(nGtype == 10);
        polymorphic = SNP_set(snp1, ndx1, 0.0345, gtype, bootchr, 1);
        assert(polymorphic);
        nGtype = encodeHaploid(gtype, sizeof(gtype), "0101010101");
        assert(nGtype == 10);
        polymorphic = SNP_set(snp2, ndx2, 0.0456, gtype, bootchr, 1);
        assert(polymorphic);
        s = gsl_rng_uniform(rng) * windowcm;
    }
    BootChr_free(boot);

    nSNPs = 100000;
    nReps = 100;
    blockLength = 300;
    nBins = 20;
    windowcm = 0.3;
    boot = BootChr_new(nSNPs, nReps, npat, blockLength,
                    windowcm, nBins, rng);
    SNP_free(snp1);
    SNP_free(snp2);
    snp1 = SNP_new(10, BootChr_nReps(boot));
    snp2 = SNP_new(10, BootChr_nReps(boot));

    /* extreme cases */
    encodeHaploid(gtype, sizeof(gtype), "0101010101");
    SNP_set(snp1, 123, 0.0345, gtype, bootchr, 1);
    encodeHaploid(gtype, sizeof(gtype), "0101010101");
    SNP_set(snp2, 234, 0.0456, gtype, bootchr, 1);
    sep = windowcm - DBL_EPSILON;
    assert(sep < windowcm);

    encodeHaploid(gtype, sizeof(gtype), "0101010101");
    SNP_set(snp1, 0, 0.0345, gtype, bootchr, 1);

    encodeHaploid(gtype, sizeof(gtype), "0101010101");
    SNP_set(snp1, Boot_nSNPs(boot) - 1, 0.0345, gtype, bootchr, 1);
    Boot_addLD(bootchr, 1.0, 1.0, s, snp1, snp2);  /*ndx1=max */
    unitTstResult("Boot_addLD", "OK");

    BootConf   *bc = BootConf_new(bootchr, 0.9);

    assert(bc);

    Boot       *boot2 = Boot_dup(boot);

    assert(Boot_equals(bootchr, boot2));
    BootConf   *bc2 = BootConf_new(boot2, 0.9);

    assert(bc2);

    Boot_free(boot2);
    boot2 = NULL;

    FILE       *dumpfile;

    /* test dump and restore */
    const char *fname = "xboot.txt";

    dumpfile = fopen(fname, "w");
    assert(dumpfile);
    Boot_dump(bootchr, dumpfile);
    fclose(dumpfile);
    dumpfile = fopen(fname, "r");
    assert(dumpfile);
    boot2 = Boot_restore(dumpfile);
    assert(Boot_equals(bootchr, boot2));
    fclose(dumpfile);
    dumpfile = NULL;
    remove(fname);

    if(verbose) {
        printf("\nthis was dumped using Boot_dump:\n");
        Boot_print(bootchr, stdout);
        BootConf_print(bc, stdout);

        printf("\nthis was restored using Boot_restore:\n");
        Boot_print(boot2, stdout);
        BootConf_print(bc2, stdout);
    }
    unitTstResult("Boot_dump", "OK");
    unitTstResult("Boot_restore", "OK");

    Boot_free(boot2);
    boot2 = NULL;

    Boot_free(boot);
    BootConf_free(bc);
    BootConf_free(bc2);
    gsl_rng_free(rng);
    SNP_free(snp1);
    SNP_free(snp2);

    unitTstResult("Boot", "OK");
#endif
    return 0;
}
