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

    int         verbose = 0;

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

    VCFReader *r = VCFReader_new("dummy.vcf");
    VCFReader_print(r, stdout);

    putchar('\n');
    VCFReader_next(r);
    VCFReader_print(r, stdout);
    VCFReader_free(r);

    unitTstResult("VCFReader", "untested");
    return 0;
}
