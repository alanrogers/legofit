/**
 * @file collapse.c
 * @author Alan R. Rogers
 * @brief Remove or merge populations from LblNdx and BranchTab.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "collapse.h"
#include "branchtab.h"
#include <stdio.h>
#include <stdlib.h>

LblNdxBranchTab removePops(LblNdx lndx, BranchTab *bt, const char *deleteStr) {
    LblNdxBranchTab rval = {
                            .lndx = lndx,
                            .branchtab = NULL
    };

    if(deleteStr == NULL) {
        fprintf(stderr,"%s:%d: NULL deleteStr\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    // Parse deleteStr to create an integer, whose "on" bits
    // indicate the populations to be deleted.
    tipId_t bits = LblNdx_getTipId(&lndx, deleteStr);
    if(bits == 0) {
        fprintf(stderr, "%s:%s:%d: delete string (%s) includes"
                " unknown labels.\n",
                __FILE__, __func__, __LINE__, deleteStr);
        fprintf(stderr,"Known labels:\n");
        LblNdx_print(&lndx, stderr);
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

LblNdxBranchTab collapsePops(LblNdx lndx, BranchTab *bt,
                             const char *collapseStr,
                             const char *lbl) {
    LblNdxBranchTab rval = {
                            .lndx = lndx,
                            .branchtab = NULL
    };

    if(collapseStr == NULL || lbl==NULL) {
        fprintf(stderr,"%s:%d: NULL deleteStr\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    // Parse collapsedStr to create an integer, whose "on" bits
    // indicate the populations to be collapsed into a single
    // population.
    tipId_t bits = LblNdx_getTipId(&lndx, collapseStr);
    if(bits == 0) {
        fprintf(stderr, "%s:%s%d: collapse string (%s) includes"
                " unknown labels.\n",
                __FILE__, __func__, __LINE__, collapseStr);
        fprintf(stderr,"Known labels:\n");
        LblNdx_print(&lndx, stderr);
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
