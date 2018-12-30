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
#include "error.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

const char *badInput =
    "#chr\tpos\tref\talt\traf\n"
    "10\t1\ta\tt\t5e-1\n"
    "1\t1\ta\t.\t0\n"
    "10\t200\tg\tc\t1e0\n"
    "10\t201\tg\tc\t1e99999\n"
    "10\t202\tg\tc\t1e-99999\n" "10\t203\tg\tc\tnot-a-float\n";

const char *tstInput[3] = {
    "#chr\tpos\tref\talt\traf\n"
    "1\t1\ta\tt\t0\n"
    "10\t1\ta\tt\t5e-1\n"
    "10\t200\tg\tc\t1e0\n",

    "#chr\tpos\tref\talt\traf\n"
    "1\t1\ta\tt\t0.5\n"
    "1\t2\ta\t.\t0.5\n"
    "10\t1\ta\tt\t1e-1\n"
    "10\t200\tg\tc\t1\n"
    "10\t201\tg\tc\t1\n",

    "#chr\tpos\tref\talt\traf\n"
    "1\t1\ta\t.\t1\n"
    "10\t1\ta\tt\t0\n"
    "10\t100\ta\tt\t0\n"
    "10\t200\tg\tc\t0\n"
};

int main(int argc, char **argv) {

    int         i, verbose = 0, status;
    char        errbuff[200];

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

    const char *badraf = "bad.raf";
    FILE       *badfp = fopen(badraf, "w");
    assert(badfp);
    fputs(badInput, badfp);
    fclose(badfp);
    RAFReader  *rdr = RAFReader_new(badraf);

    for(i=0; i<7; ++i) {
        status = RAFReader_next(rdr);
        if(i==0)
            assert(status==0);
        else if(i==1)
            assert(status==BAD_SORT);
        else if(i==2)
            assert(status==0);
        else if(i==3)
            assert(status==ERANGE);
        else if(i==4)
            assert(status==ERANGE);
        else if(i==5)
            assert(status==EINVAL);
        else if(i==6)
            assert(status==EOF);
    }while(status != EOF);
    RAFReader_free(rdr);
    remove(badraf);

    const char *tst[3] = {"tst0.raf", "tst1.raf", "tst2.raf"};
    FILE *      fp[3];
    for(i = 0; i < 3; ++i) {
        fp[i] = fopen(tst[i], "w");
        assert(fp[i]);
        fputs(tstInput[i], fp[i]);
        fclose(fp[i]);
    }

    RAFReader * r[3];
    if(verbose)
        RAFReader_printHdr(stderr);
    for(i = 0; i < 3; ++i) {
        r[i] = RAFReader_new(tst[i]);
        do{
            status = RAFReader_next(r[i]);
            switch (status) {
            case 0:
                if(verbose)
                    RAFReader_print(r[i], stderr);
                break;
            case EOF:
                break;
            default:
                mystrerror_r(status, errbuff, sizeof errbuff);
                fprintf(stderr, "%s:%d: i=%d %s\n",
                        __FILE__,__LINE__, i, errbuff);
                exit(1);
            }
        }while(status != EOF);
        RAFReader_free(r[i]);
    }

    for(i = 0; i < 3; ++i)
        r[i] = RAFReader_new(tst[i]);

    i=0;
    do{
        status = RAFReader_multiNext(3, r);
        if(status==0)
            status = RAFReader_findDaf(3, r);
        if(i==0) {
            assert(status==0);
            assert(0 == strcmp("1",RAFReader_chr(r[0])));
            assert(1UL == RAFReader_nucpos(r[0]));
            assert(1.0 == RAFReader_daf(r[0]));
            assert(0.5 == RAFReader_daf(r[1]));
            assert(0.0 == RAFReader_daf(r[2]));
        }else if(i==1) {
            assert(status==0);
            assert(0 == strcmp("10",RAFReader_chr(r[0])));
            assert(1UL == RAFReader_nucpos(r[0]));
            assert(0.5 == RAFReader_daf(r[0]));
            assert(0.1 == RAFReader_daf(r[1]));
            assert(0.0 == RAFReader_daf(r[2]));
        }else if(i==2) {
            assert(status==NO_ANCESTRAL_ALLELE);
            assert(0 == strcmp("10",RAFReader_chr(r[0])));
            assert(200UL == RAFReader_nucpos(r[0]));
        }else if(i==3) {
            assert(status==EOF);
        }else{
            fprintf(stderr,"%s:%d: this shouldn't happen\n",
                    __FILE__,__LINE__);
            exit(1);
        }
        ++i;
    }while(status != EOF);

    for(i = 0; i < 3; ++i) {
        RAFReader_free(r[i]);
        //remove(tst[i]);
    }

    unitTstResult("RAFReader", "OK");

    return 0;
}
