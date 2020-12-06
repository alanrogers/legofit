#include "collapse.h"
#include "branchtab.h"
#include <stdio.h>
#include <stdlib.h>

LblNdxBranchTab removePops(LblNdx lndx, BranchTab *bt, char *deleteStr) {
    LblNdxBranchTab rval;

    // Parse deleteStr to create an integer, whose "on" bits
    // indicate the populations to be deleted.
    tipId_t bits = LblNdx_getTipId(&lndx, deleteStr);
    if(bits == 0) {
        fprintf(stderr, "%s:%d: delete string (%s) includes"
                " nonexistent populations.\n",
                __FILE__, __LINE__, deleteStr);
        exit(EXIT_FAILURE);
    }

    // Create a new LblNdx object within rval, which lacks the
    // populations we are removing.
    int status = LblNdx_rmPops(&(rval.lndx), bits);
    if(status) {
        fprintf(stderr,"%s:%d: can't remove a population that"
                " doesn't exist.\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    // Create a new BranchTab w/i rval, which lack deleted populations.
    rval.branchtab = BranchTab_rmPops(bt, bits);

    return rval;
}

LblNdxBranchTab collapsePops(LblNdx lndx, BranchTab *bt, char *collapseStr,
                         const char *lbl) {
    LblNdxBranchTab rval;

    // Parse collapsedStr to create an integer, whose "on" bits
    // indicate the populations to be collapsed into a single
    // population.
    tipId_t bits = LblNdx_getTipId(&lndx, collapseStr);
    if(bits == 0) {
        fprintf(stderr, "%s:%d: collapse string (%s) includes"
                " nonexistent populations.\n",
                __FILE__, __LINE__, collapseStr);
        exit(EXIT_FAILURE);
    }

    // Create a new LblNdx object within rval, with several populations
    // collapsed into one.
    int status = LblNdx_collapse(&(rval.lndx), bits, lbl);
    if(status) {
        fprintf(stderr,"%s:%d: string %s has a population that"
                " doesn't exist.\n",
                __FILE__,__LINE__, collapseStr);
        exit(EXIT_FAILURE);
    }

    // Create a new BranchTab w/i rval, which lack deleted populations.
    rval.branchtab = BranchTab_collapse(bt, bits);

    return rval;
}
