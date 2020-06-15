/**
 * @file xmigoutcome.c
 * @author Alan R. Rogers
 * @brief Test migoutcome.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "migoutcome.h"
#include "misc.h"
#include <stdio.h>
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
            fprintf(stderr, "usage: xmigoutcome [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xmigoutcome [-v]\n");
        exit(1);
    }

    MigOutcome *mo = NULL;

    for(unsigned i=0; i < 3; ++i) {
        unsigned event = nextMigrationEvent();
        assert(i == event);
        unsigned noutcomes = (i+1) << i;
        for(unsigned outcome=0; outcome < noutcomes; ++outcome)
            mo = MigOutcome_insert(mo, event, noutcomes, outcome, 0.5);
    }

    MigOutcome *mo2 = MigOutcome_dup(mo);
    if(verbose) {
        MigOutcome_print(mo, stdout);
        MigOutcome_print(mo2, stdout);
    }

    MigOutcome *a=mo, *b=mo2;
    while(a && b) {
        assert(a->event == b->event);
        assert(a->noutcomes == b->noutcomes);
        assert(a->bits == b->bits);
        assert(a->pr == b->pr);
        a = a->next;
        b = b->next;
    }
    assert(a==NULL);
    assert(b==NULL);
    MigOutcome_free(mo);
    MigOutcome_free(mo2);

    unitTstResult("MigOutcome", "OK");
    
    return 0;
}
