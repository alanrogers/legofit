/**
 * @file xvcfreader.c
 * @author Alan R. Rogers
 * @brief Test vcfreader.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "vcfreader.h"
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
            fprintf(stderr, "usage: xvcfreader [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xvcfreader [-v]\n");
        exit(1);
    }

    VCFReader *r[3];
	r[0] = VCFReader_new("dummy.vcf");
	while(EOF != VCFReader_next(r[0])) {
		if(verbose)
			VCFReader_print(r[0], stdout);
	}
    VCFReader_free(r[0]);

	r[0] = VCFReader_new("dummy.vcf");
	r[1] = VCFReader_new("dummy2.vcf");
	r[2] = VCFReader_new("dummy3.vcf");

	while(EOF != VCFReader_multiNext(3, r)) {
		if(verbose) {
			fputs("#################\n", stdout);
			for(i=0; i<3; ++i)
				VCFReader_print(r[i], stdout);
		}
	}

    unitTstResult("VCFReader", "untested");
    return 0;
}
