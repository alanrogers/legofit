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

    RAFReader *r[3];
	r[0] = RAFReader_new("altai.raf");
	while(EOF != RAFReader_next(r[0])) {
		//      const char *chr = RAFReader_chr(r[0]);
        //		assert(0 == strcmp(chr, "22"));
		if(verbose)
			RAFReader_print(r[0], stdout);
	}
    RAFReader_free(r[0]);

	r[0] = RAFReader_new("altai.raf");
	r[1] = RAFReader_new("denisova.raf");
	r[2] = RAFReader_new("Mgenomes3.raf");

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
           100* mismatch / ((double) (match+mismatch)));

    RAFReader_free(r[0]);
    RAFReader_free(r[1]);
    RAFReader_free(r[2]);

    unitTstResult("RAFReader", "untested");
    return 0;
}
