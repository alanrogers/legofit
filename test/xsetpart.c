/**
 * @file xsetpart.c
 * @author Alan R. Rogers
 * @brief Generate all partitions of a set of n elements into m subsets.
 * @copyright Copyright (c) 2019, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "setpart.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

void usage(void);
int visit(unsigned n, unsigned a[n], void *data);

typedef struct VisitDat VisitDat;

// Data manipulated by visit function.
struct VisitDat {
    int verbose;
    unsigned nsubdivisions;
    unsigned long count;
    double lnconst; // log of constant in Durrett's theorem 1.5
    double sumprob; // sum of probabilities should equal 1.
};

/// Function to process partition. This function prints
/// the partition if verbose!=0, checks that it is valid,
/// and counts the partitions. Returns non-zero if an error is found.
int visit(unsigned n, unsigned a[n], void *data) {
    VisitDat *vdat = (VisitDat *) data;

    unsigned max = 0;
    int status=0;
    for(int i=0; i<n; ++i) {
        if(a[i] > max) {
            if(a[i] != max+1)
                ++status;
            max = a[i];
        }
    }
    if(max != 1 + vdat->nsubdivisions)
        ++status;
    vdat->count += 1;

    // Calculate the number, k, of partitions and the size, c[i], of each.
    unsigned k = 0;
    for(int i=0; i<n; ++i) {
        if(a[i] > k)
            k = a[i];
    }
    ++k;
    unsigned c[k];
    memset(c, 0, k*sizeof(c[0]));
    for(int i=0; i<n; ++i)
        ++c[a[i]];

    double prob = probPartition(k, c, vdat->lnconst);
    vdat->sumprob += prob;

    if(vdat->verbose) {
        for(int i=0; i<n; ++i) {
            printf("%d", a[i]);
            if(i < n-1)
                putchar(' ');
        }
        printf(" : %le\n", prob);
    }

    return status;
}

void usage(void) {
    fprintf(stderr,"usage: xsetpart [options]\n");
    fprintf(stderr,"where options may include\n");
    fprintf(stderr," -n <n_elements> : size of set\n");
    fprintf(stderr," -k <n_subdiv>   : number of subdivisions\n");
    fprintf(stderr," -s              : test Stirling2 only\n");
    fprintf(stderr," -v              : verbose output\n\n");
    fprintf(stderr,"   n_elements must be > 0\n");
    fprintf(stderr,"   n_subdiv must be > 0, <= n_elements\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
    int i, n, verbose = 0, s2only=0;
    unsigned nelem = 6, nsub = 3;

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose=1;
        else if(strcmp(argv[i], "-n") == 0) {
            i += 1;
            if(i == argc)
                usage();
            nelem = strtoul(argv[i], NULL, 10);
        }else if(strcmp(argv[i], "-k") == 0) {
            i += 1;
            if(i == argc)
                usage();
            nsub = strtoul(argv[i], NULL, 10);
        }else if(strcmp(argv[i], "-s") == 0) {
            s2only = 1;
        }else
            usage();
    }

    Stirling2 *s = Stirling2_new(nelem);
    CHECKMEM(s);

    if(verbose) {
        Stirling2_print(s, stdout);
        long unsigned z = Stirling2_val(s, nelem, nsub);
        printf("Stirling2(%u, %u) = %lu = %le\n",
               nelem, nsub, z, (double) z);
    }

    if(s2only)
        return 0;

    assert(Stirling2_val(s, 0, 0) == 1);
    for(n=1; n <= nelem; ++n) {
        assert(Stirling2_val(s,n,0) == 0);
        assert(Stirling2_val(s,n,1) == 1);
        assert(Stirling2_val(s,n,n) == 1);
    }

    if(nelem >= 3)
        assert(Stirling2_val(s, 3, 2) == 3);

    if(nelem >= 4) {
        assert(Stirling2_val(s, 4, 2) == 7);
        assert(Stirling2_val(s, 4, 3) == 6);
    }

    if(nelem >= 5) {
        assert(Stirling2_val(s, 5, 2) == 15);
        assert(Stirling2_val(s, 5, 3) == 25);
        assert(Stirling2_val(s, 5, 4) == 10);
    }

    if(nelem >= 6) {
        assert(Stirling2_val(s, 6, 2) == 31);
        assert(Stirling2_val(s, 6, 3) == 90);
        assert(Stirling2_val(s, 6, 4) == 65);
        assert(Stirling2_val(s, 6, 5) == 15);
    }

    Stirling2_free(s);

    unitTstResult("Stirling2", "OK");

    VisitDat vdat = {.verbose = verbose, 
                     .nsubdivisions = nsub,
                     .count=0,
                     .sumprob=0.0,
                     .lnconst=lnCoalConst(nelem, nsub) };

    if(verbose)
        printf("lnconst=%lf\n", vdat.lnconst);

    int status = traverseSetPartitions(nelem, nsub, visit, &vdat);

    if(verbose)
        printf("sumprob=%lf; should be 1.0; err=%le\n",
               vdat.sumprob, vdat.sumprob - 1.0);
    assert(fabs(vdat.sumprob - 1.0) < 1e-8);

    s = Stirling2_new(nelem);
    CHECKMEM(s);
    if(verbose)
        printf("Stirling2(%u,%u)=%lu\n", nelem, nsub,
               Stirling2_val(s, nelem, nsub));
    assert(vdat.count == Stirling2_val(s, nelem, nsub));

    unitTstResult("SetPart", status==0 ? "OK" : "FAIL");
    return 0;
}
