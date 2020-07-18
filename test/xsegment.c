/**
 * @file xsegment.c
 * @author Alan R. Rogers
 * @brief Test segment.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "segment.h"
#include "misc.h"
#include "error.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char *argv[]) {

    int verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            goto usage;
        verbose = 1;
        break;
    default:
        goto usage;
    }

    int         i, status, ok=1;
    const int   nseg = 7;
    Segment     v[nseg];
    Segment    *s[nseg];
    NodeStore  *ns = NodeStore_new(nseg, sizeof(v[0]), v);
    double      twoN=100.0;
    double      t0=0.0, t1=1.0, t2=2.0, t3=3.0;
    double      mix = 0.5;

    CHECKMEM(ns);


    s[0] = Segment_new(&twoN, &t0, 0, ns);
    s[1] = Segment_new(&twoN, &t0, 0, ns);
    s[2] = Segment_new(&twoN, &t0, 0, ns);
    s[3] = Segment_new(&twoN, &t1, 0, ns);
    s[4] = Segment_new(&twoN, &t1, 0, ns);
    s[5] = Segment_new(&twoN, &t2, 0, ns);
    s[6] = Segment_new(&twoN, &t3, 0, ns);

    for(i=0; i<nseg; ++i) {
        assert(s[i]->twoN == &twoN);
        assert(s[i]->end == NULL);
        assert(s[i]->mix == NULL);
        assert(s[i]->nsamples == 0);
        assert(s[i]->nchildren == 0);
        assert(s[i]->child[0] == NULL);
        assert(s[i]->child[1] == NULL);
        assert(s[i]->parent[0] == NULL);
        assert(s[i]->parent[1] == NULL);
    }

    NodeStore_free(ns);

    char buff[100];

    /*
      2-4----+
        |    |
        V    6-
      1-3-+  |
          5--+
      0---+
     */

    // 1 is a mixture of 3 and 4
    status = Segment_mix(s[1], &mix, s[4], s[3]);
    if(status) {
        mystrerror_r(status, buff, sizeof(buff));
        fprintf(stderr,"%s:%d: %s\n", __FILE__,__LINE__,buff);
        ok=0;
    }

    status = Segment_addChild(s[5], s[0]); // derive 0 from 5
    if(status) {
        mystrerror_r(status, buff, sizeof(buff));
        fprintf(stderr,"%s:%d: %s\n", __FILE__,__LINE__,buff);
        ok=0;
    }

    status = Segment_addChild(s[5], s[3]); // derive 3 from 5
    if(status) {
        mystrerror_r(status, buff, sizeof(buff));
        fprintf(stderr,"%s:%d: %s\n", __FILE__,__LINE__,buff);
        ok=0;
    }

    status = Segment_addChild(s[4], s[2]); // derive 2 from 4
    if(status) {
        mystrerror_r(status, buff, sizeof(buff));
        fprintf(stderr,"%s:%d: %s\n", __FILE__,__LINE__,buff);
        ok=0;
    }

    status = Segment_addChild(s[6], s[4]); // derive 4 from 6
    if(status) {
        mystrerror_r(status, buff, sizeof(buff));
        fprintf(stderr,"%s:%d: %s\n", __FILE__,__LINE__,buff);
        ok=0;
    }

    status = Segment_addChild(s[6], s[5]); // derive 5 from 6
    if(status) {
        mystrerror_r(status, buff, sizeof(buff));
        fprintf(stderr,"%s:%d: %s\n", __FILE__,__LINE__,buff);
        ok=0;
    }

    Segment *root = Segment_root(s[0]);
    for(i=1; i<nseg; ++i) {
        Segment *r = Segment_root(s[i]);
        if(r != root) {
            fprintf(stderr,"%s:%d segments 0 and %d have different roots\n",
                    __FILE__,__LINE__, i);
            ok=0;
        }
    }

    Segment_sanityCheck(root, __FILE__, __LINE__);

    if(verbose)
        Segment_print(stdout, root, 0);
    
    unitTstResult("Segment", ok ? "OK": "FAIL");

    return 0;

 usage:
    fprintf(stderr, "usage xsegment [-v]\n");
    exit(EXIT_FAILURE);
}
