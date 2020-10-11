/**
 * @file xeventlst.c
 * @author Alan R. Rogers
 * @brief Test event.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "eventlst.h"
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
            fprintf(stderr, "usage: xevent [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xevent [-v]\n");
        exit(1);
    }

    EventLst *mo[2] = {NULL, NULL};

    for(unsigned i=0; i < 3; ++i) {
        unsigned event = nextEvent();
        for(unsigned outcome = 0; outcome < 2; ++outcome) {
            mo[outcome] = EventLst_insert(mo[outcome], event,
                                            outcome, 0.5);
        }
    }

    assert(0 > EventLst_cmp(mo[0], mo[1]));

    EventLst *mo3 = EventLst_dup(mo[0]);
    if(verbose) {
        EventLst_print(mo[0], stdout);
        EventLst_print(mo[1], stdout);
        EventLst_print(mo3, stdout);
    }

    assert(0 == EventLst_cmp(mo[0], mo3));

    int mutually_exclusive;
    EventLst *mo4 = EventLst_join(mo[0], mo[1], &mutually_exclusive);
    assert(mo4 == NULL);
    assert(mutually_exclusive == 1);

    EventLst_free(mo3);
    mo3 = NULL;
    for(unsigned i=0; i < 2; ++i) {
        unsigned event = nextEvent();
        mo3 = EventLst_insert(mo3, event, 1, 0.5);
    }
    mo4 = EventLst_join(mo[0], mo3, &mutually_exclusive);
    assert(mo4);
    assert(mutually_exclusive == 0);
    if(verbose) {
        EventLst_print(mo4, stdout);
        putchar('\n');
    }

    EventLst_free(mo4);
    mo4 = NULL;
    mo4 = EventLst_join(mo[0], mo[0], &mutually_exclusive);
    assert(mo4);
    assert(mutually_exclusive == 0);

    EventLst *a=mo[0], *b=mo[1];
    while(a && b) {
        assert(a->event == b->event);
        assert(a->outcome != b->outcome);
        assert(a->pr == b->pr);
        a = a->next;
        b = b->next;
    }
    assert(a==NULL);
    assert(b==NULL);
    EventLst_free(mo[0]);
    EventLst_free(mo[1]);
    EventLst_free(mo3);
    EventLst_free(mo4);

    unitTstResult("EventLst", "OK");
    
    return 0;
}
