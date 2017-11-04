/**
 * @file xrafreader.c
 * @author Alan R. Rogers
 * @brief Test rafreader.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "rafreader.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

const char *tstInput[3] = {"#chr\tpos\tref\talt\traf\n"
                           "1\t1\ta\t.\t0\n"
                           "10\t1\ta\tt\t5e-1\n"
                           "10\t200\tg\tc\t1e0\n",

                           "#chr\tpos\tref\talt\traf\n"
                           "1\t1\ta\t.\t0.5\n"
                           "1\t2\ta\t.\t0.5\n"
                           "10\t1\ta\tt\t1e-1\n"
                           "10\t200\tg\tc\t1\n",

                           "#chr\tpos\tref\talt\traf\n"
                           "1\t1\ta\t.\t0.123\n"
                           "10\t1\ta\tt\t5e-1\n"
                           "10\t100\ta\tt\t5e-3\n"
                           "10\t200\tg\tc\t0.000\n"};

int main(int argc, char **argv) {

    int         i, verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xrafreader [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xrafreader [-v]\n");
        exit(1);
    }

    const char *tst[3] = {"tst0.raf", "tst1.raf", "tst2.raf"};
    FILE *fp[3];
    for(i=0; i<3; ++i) {
        fp[i] = fopen(tst[i], "w");
        assert(fp[i]);
        fputs(tstInput[i], fp[i]);
        fclose(fp[i]);
    }

    RAFReader *r[3];
    for(i=0; i<3; ++i) {
        r[i] = RAFReader_new(tst[i]);
        while(EOF != RAFReader_next(r[i])) {
            if(verbose)
                RAFReader_print(r[i], stdout);
        }
        RAFReader_free(r[0]);
    }

    for(i=0; i<3; ++i)
        r[i] = RAFReader_new(tst[i]);

    long match=0, mismatch=0;
	while(EOF != RAFReader_multiNext(3, r)) {
		if(verbose) {
			fputs("#################\n", stdout);
			for(i=0; i<3; ++i)
				printf(" raf[%d]=%lf", i, RAFReader_raf(r[i]));
            putchar('\n');
		}
        if(RAFReader_allelesMatch(3, r))
            ++match;
        else {
            printf("Mismatch\n");
            RAFReader_printHdr(stdout);
            for(i=0; i<3; ++i)
                RAFReader_print(r[i], stdout);
            ++mismatch;
        }
	}
    printf("%ld/%ld (%lf%%) of SNPs have ref/alt alleles that don't match.\n",
           mismatch, match+mismatch,
           100*mismatch / ((double) (match+mismatch)));

    for(i=0; i<3; ++i) {
        RAFReader_free(r[i]);
        remove(tst[i]);
    }

    unitTstResult("RAFReader", "untested");

    return 0;
}
