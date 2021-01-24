/**
 * @file xstrdblqueue.c
 * @author Daniel R. Tabin
 * @brief Unit tests for StrDblQueue
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "hessian.h"
#include "misc.h"
#include "strdblqueue.h"
#include <stdbool.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv){
    StrDblQueue* a = NULL;
    StrDblQueue* b = NULL;
    StrDblQueue* c = NULL;
    StrDblQueue* d = NULL;
    StrDblQueue* e = NULL;
    StrDblQueue* f = NULL;

    const char *fakeBoot =
    "########################################\n"
    "# legofit: estimate population history #\n"
    "########################################\n"
    "\n"
    "# Program was compiled: Oct  7 2017 11:29:11\n"
    "# Program was run: Sat Oct  7 15:09:02 2017\n"
    ""
    "# cmd: legofit -1 --tol 5e-4 -s 4 -S 1000@10000 -S 2000@2000000 s1.lgo ../boot/sing048\n"
    "# Stage nOptItr nSimReps\n"
    "#     0    1000    10000\n"
    "#     1    2000  2000000\n"
    "# DE strategy        : 4\n"
    "#    F               : 0.900000\n"
    "#    CR              : 0.800000\n"
    "#    tolerance       : 5.000000e-04\n"
    "# nthreads           : 42\n"
    "# lgo input file     : s1.lgo\n"
    "# site pat input file: ../boot/sing048\n"
    "# pts/dimension      : 10\n"
    "# Including singleton site patterns.\n"
    "# cost function      : KL\n"
    "Initial parameter values\n"
    "    4 fixed:\n"
    "       zero = 0\n"
    "        one = 1\n"
    "      Txynd = 25920\n"
    "        TmN = 1897\n"
    "    0 Gaussian:\n"
    "    9 free:\n"
    "        Tnd = 17318\n"
    "        Txy = 3788\n"
    "         Td = 1734\n"
    "         Ta = 1760\n"
    "       2Nnd = 1000\n"
    "        2Nn = 1000\n"
    "       2Nxy = 36808\n"
    "     2Nxynd = 35097\n"
    "         mN = 0.0246\n"
    "    0 constrained:\n"
    "DiffEv converged. cost=5.30424e-04 spread=4.96725e-04\n"
    "\n"
    "Fitted parameter values\n"
    "    9 free:\n"
    "        Tnd = 21907.8\n"
    "        Txy = 9836.87\n"
    "         Td = 782.991\n"
    "         Ta = 1852.97\n"
    "       2Nnd = 9994.27\n"
    "        2Nn = 20166.8\n"
    "       2Nxy = 31661.7\n"
    "     2Nxynd = 45033.9\n"
    "         mN = 0.030938\n"
    "    0 constrained:\n"
    "#       SitePat  BranchLen\n"
    "              x 38177.7670041\n"
    "              y 36976.1129329\n"
    "              n 39773.2585902\n"
    "              d 41503.9037230\n"
    "            x:y 22741.5617898\n"
    "            x:n 2973.2394928\n"
    "            x:d 3177.9284049\n"
    "            y:n 3699.1331968\n"
    "            y:d 3002.9907387\n"
    "            n:d 19125.0587195\n"
    "          x:y:n 7062.0871479\n"
    "          x:y:d 6892.8717880\n"
    "          x:n:d 5889.0132020\n"
    "          y:n:d 6539.7112355\n";

    FILE       *fakeFile = fopen("s1boot0.legofit", "w");
    assert(fakeFile);
    fputs(fakeBoot, fakeFile);
    fclose(fakeFile);

    StrDbl temp;

    bool verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xstrdblqueue [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xstrdblqueue [-v]\n");
        exit(EXIT_FAILURE);
    }

    assert(StrDblQueue_compare(a,b) == 0);
    assert(StrDblQueue_length(a) == 0);

    a = StrDblQueue_push(a, "test", 1.0);

    if(verbose){
      fprintf(stderr,"Queue a:\n");
      StrDblQueue_print(a,stderr);
      fprintf(stderr,"Queue b:\n");
      StrDblQueue_print(b,stderr);
    }

    assert(StrDblQueue_compare(a,b) != 0);
    assert(StrDblQueue_length(a) == 1);

    a = StrDblQueue_pop(a, &temp);
    a = StrDblQueue_pop(a, &temp);

    if(verbose){
      fprintf(stderr,"Queue a:\n");
      StrDblQueue_print(a,stderr);
      fprintf(stderr,"Queue b:\n");
      StrDblQueue_print(b,stderr);
    }

    assert(StrDblQueue_compare(a,b) == 0);
    assert(StrDblQueue_length(a) == 0);

    c = StrDblQueue_parseLegofit("s1boot0.legofit");
    d = StrDblQueue_parseLegofit("s1boot0.legofit");

    assert(StrDblQueue_compare(c,d) == 0);

    if(verbose){
      fprintf(stderr,"Queue c:\n");
      StrDblQueue_print(c,stderr);
      fprintf(stderr,"Queue d:\n");
      StrDblQueue_print(d,stderr);
    }

    c = StrDblQueue_push(c, "test", 1);

    if(verbose){
      fprintf(stderr,"Queue c:\n");
      StrDblQueue_print(c,stderr);
      fprintf(stderr,"Queue d:\n");
      StrDblQueue_print(d,stderr);
    }

    assert(StrDblQueue_compare(c,d) != 0);

    c = StrDblQueue_pop(c, &temp);
    d = StrDblQueue_pop(d, &temp);

    assert(StrDblQueue_compare(c,d) != 0);

    if(verbose){
      fprintf(stderr,"Queue c:\n");
      StrDblQueue_print(c,stderr);
      fprintf(stderr,"Queue d:\n");
      StrDblQueue_print(d,stderr);
    }

    e = StrDblQueue_parseSitePat("s1boot0.legofit");
    f = StrDblQueue_parseSitePat("s1boot0.legofit");

    assert(StrDblQueue_compare(e,f) == 0);

    if(verbose){
      fprintf(stderr,"Queue e:\n");
      StrDblQueue_print(e,stderr);
      fprintf(stderr,"Queue f:\n");
      StrDblQueue_print(f,stderr);
    }

    e = StrDblQueue_push(e, "test", 1);

    if(verbose){
      fprintf(stderr,"Queue e:\n");
      StrDblQueue_print(e,stderr);
      fprintf(stderr,"Queue f:\n");
      StrDblQueue_print(f,stderr);
    }

    assert(StrDblQueue_compare(e,f) != 0);

    e = StrDblQueue_pop(e, &temp);
    f = StrDblQueue_pop(f, &temp);

    assert(StrDblQueue_compare(e,f) != 0);

    if(verbose){
      fprintf(stderr,"Queue e:\n");
      StrDblQueue_print(e,stderr);
      fprintf(stderr,"Queue f:\n");
      StrDblQueue_print(f,stderr);
    }

    remove("s1boot0.legofit");

    a = StrDblQueue_free(a);
    b = StrDblQueue_free(b);
    c = StrDblQueue_free(c);
    d = StrDblQueue_free(d);
    e = StrDblQueue_free(e);
    f = StrDblQueue_free(f);

    assert(a==NULL);
    assert(b==NULL);
    assert(c==NULL);
    assert(d==NULL);
    assert(e==NULL);
    assert(f==NULL);

    a = StrDblQueue_push(a, "x", 1.0);
    b = StrDblQueue_push(b, "x", 2.0);
    a = StrDblQueue_push(a, "y", 3.0);
    b = StrDblQueue_push(b, "y", 4.0);
    assert(1.0 == StrDblQueue_msd(a,b));

    a = StrDblQueue_free(a);
    b = StrDblQueue_free(b);
    c = StrDblQueue_free(c);
    d = StrDblQueue_free(d);
    e = StrDblQueue_free(e);
    f = StrDblQueue_free(f);

    unitTstResult("StrDblQueue", "OK");

    return 0;
}
