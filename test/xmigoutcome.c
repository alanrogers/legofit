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

    MigOutcome *mo[2] = {NULL, NULL};

    for(unsigned i=0; i < 3; ++i) {
        unsigned event = nextMigrationEvent();
        for(unsigned outcome = 0; outcome < 2; ++outcome) {
            mo[outcome] = MigOutcome_insert(mo[outcome], event,
                                            outcome, 0.5);
        }
    }

    MigOutcome *mo3 = MigOutcome_dup(mo[0]);
    if(verbose) {
        MigOutcome_print(mo[0], stdout);
        MigOutcome_print(mo[1], stdout);
        MigOutcome_print(mo3, stdout);
    }

    int mutually_exclusive;
    MigOutcome *mo4 = MigOutcome_join(mo[0], mo[1], &mutually_exclusive);
    assert(mo4 == NULL);
    assert(mutually_exclusive == 1);

    MigOutcome_free(mo3);
    mo3 = NULL;
    for(unsigned i=0; i < 2; ++i) {
        unsigned event = nextMigrationEvent();
        mo3 = MigOutcome_insert(mo3, event, 1, 0.5);
    }
    mo4 = MigOutcome_join(mo[0], mo3, &mutually_exclusive);
    assert(mo4);
    assert(mutually_exclusive == 0);
    if(verbose)
        MigOutcome_print(mo4, stdout);

    MigOutcome_free(mo4);
    mo4 = NULL;
    mo4 = MigOutcome_join(mo[0], mo[0], &mutually_exclusive);
    assert(mo4);
    assert(mutually_exclusive == 0);

    MigOutcome *a=mo[0], *b=mo[1];
    while(a && b) {
        assert(a->event == b->event);
        assert(a->outcome != b->outcome);
        assert(a->pr == b->pr);
        a = a->next;
        b = b->next;
    }
    assert(a==NULL);
    assert(b==NULL);
    MigOutcome_free(mo[0]);
    MigOutcome_free(mo[1]);
    MigOutcome_free(mo3);
    MigOutcome_free(mo4);

    unitTstResult("MigOutcome", "OK");
    
    return 0;
}
