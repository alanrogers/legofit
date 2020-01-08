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
double lnCoalConst(unsigned n, unsigned k);
double prPartition(unsigned k, unsigned y[k], double lnconst);

typedef struct VisitDat VisitDat;

// Data manipulated by visit function.
struct VisitDat {
    int verbose;
    unsigned nsubdivisions;
    unsigned long count;
    double lnconst; // log of constant in Durrett's theorem 1.5
    double sumprob; // sum of probabilities should equal 1.
};

/** 
 * Log of constant in coalescent probability from theorem 1.5, p. 11, Durrett,
 * Richard. 2008. Probability Models for DNA Sequence Evolution.
 */
double lnCoalConst(unsigned n, unsigned k) {
    assert(n >= k);
    assert(k > 0);
    double ans = lgamma(k+1)
        - lgamma(n+1)
        + lgamma(n-k+1)
        + lgamma(k-1)
        - lgamma(n-1);
    return ans;
}

/**
 * Partition probability under the coalescent process. There are n
 * descendants in some recent epoch and k < n in some earlier
 * epoch. In that earlier epoch, the i'th ancestor had y[i]
 * descendants. The function returns the probability of a partition
 * that satisfies this condition, given n and k. There may be several
 * such partitions. lnconst should be calculated using function
 * lnCoalConst. See theorem 1.5, p. 11, Durrett,
 * Richard. 2008. Probability Models for DNA Sequence Evolution.
 */
double prPartition(unsigned k, unsigned y[k], double lnconst) {
    assert(k > 0);
    double ans = lnconst;
    for(unsigned i=0; i<k; ++i)
        ans += lgamma(y[i] + 1);
    return exp(ans);
}

/// Function to process partition. This function prints
/// the partition if verbose!=0, checks that it is valid,
/// and counts the partitions. Returns non-zero if an error is found.
int visit(unsigned n, unsigned a[n], void *data) {
    VisitDat *vdat = (VisitDat *) data;
    if(vdat->verbose) {
        for(int i=0; i<n; ++i) {
            printf("%d", a[i]);
            if(i < n-1)
                putchar(' ');
        }
        putchar('\n');
    }
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
    for(int i=1; i<n; ++i) {
        ++c[a[i]];
    }
    
    vdat->sumprob += prPartition(k, c, vdat->lnconst);
    return status;
}

void usage(void) {
    fprintf(stderr,"usage: xsetpart [options]\n");
    fprintf(stderr,"where options may include\n");
    fprintf(stderr," -n <n_elements> : size of set\n");
    fprintf(stderr," -k <n_subdiv>   : number of subdivisions\n");
    fprintf(stderr," -v              : verbose output\n\n");
    fprintf(stderr,"   n_elements must be > 0\n");
    fprintf(stderr,"   n_subdiv must be > 0, <= n_elements\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
    int i, n, verbose = 0;
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
        }else
            usage();
    }

    Stirling2 *s = Stirling2_new(nelem);
    CHECKMEM(s);

    if(verbose)
        Stirling2_print(s, stdout);

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
        printf("sumprob=%lf; should by 1.0; err=%le\n",
               vdat.sumprob, vdat.sumprob - 1.0);
    assert(fabs(vdat.sumprob - 1.0) < 1e-8);

    s = Stirling2_new(nelem);
    CHECKMEM(s);
    if(verbose)
        printf("Stirling3(%u,%u)=%lu\n", nelem, nsub,
               Stirling2_val(s, nelem, nsub));
    assert(vdat.count == Stirling2_val(s, nelem, nsub));

    unitTstResult("SetPart", status==0 ? "OK" : "FAIL");
    return 0;
}
