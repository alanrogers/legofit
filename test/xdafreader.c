/**
 * @file xdafreader.c
 * @author Alan R. Rogers
 * @brief Test dafreader.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "dafreader.h"
#include "misc.h"
#include "strint.h"
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
            fprintf(stderr, "usage: xdafreader [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xdafreader [-v]\n");
        exit(1);
    }

    StrInt *strint = StrInt_new();
	StrInt_insert(strint, "22", 0);

    DAFReader *r[3];
	r[0] = DAFReader_new("altai.daf");
	while(EOF != DAFReader_next(r[0])) {
		const char *chr = DAFReader_chr(r[0]);
		assert(0 == strcmp(chr, "22"));
		assert(0 == DAFReader_chrNdx(r[0], strint));
		if(verbose)
			DAFReader_print(r[0], stdout);
	}
	StrInt_free(strint);
    DAFReader_free(r[0]);

	r[0] = DAFReader_new("altai.daf");
	r[1] = DAFReader_new("denisova.daf");
	r[2] = DAFReader_new("Mgenomes3.daf");

    strint = StrInt_new();
    for(i=0; i < 22; ++i) {
        char foo[10];
        snprintf(foo, sizeof foo, "%d", i+1);
        StrInt_insert(strint, foo, i);
    }

    long match=0, mismatch=0;
	while(EOF != DAFReader_multiNext(3, r, strint)) {
		if(verbose) {
			fputs("#################\n", stdout);
			for(i=0; i<3; ++i)
				printf(" daf[%d]=%lf", i, DAFReader_daf(r[i]));
            putchar('\n');
		}
        if(DAFReader_allelesMatch(3, r))
            ++match;
        else {
            printf("Mismatch\n");
            DAFReader_printHdr(stdout);
            for(i=0; i<3; ++i)
                DAFReader_print(r[i], stdout);
            ++mismatch;
        }
	}
    printf("%ld/%ld (%lf%%) of SNPs have ancestral/derived"
           " alleles that don't match.\n",
           mismatch, match+mismatch,
           100* mismatch / ((double) (match+mismatch)));

    StrInt_free(strint);

    unitTstResult("DAFReader", "untested");
    return 0;
}
