#include "pointbuff.h"
#include "misc.h"
#include <string.h>
#include <stdio.h>
int main(int argc, char **argv) {

    const unsigned npar = 2;
    unsigned totpts = 2;
    PointBuff *pb = PointBuff_new(npar, totpts);
    if(pb == NULL) {
        fprintf(stderr,"%s:%d: bad return from PointBuff\n",
                __FILE__,__LINE__);
        exit(1);
    }

    double cx=100.0, cy=200.0, cz = 300.0;
    double x[npar], y[npar], z[npar];
    double initval=0.0;
    int i;
    for(i=0; i<npar; ++i)
        x[i] = ++initval;
    for(i=0; i<npar; ++i)
        y[i] = ++initval;
    for(i=0; i<npar; ++i)
        z[i] = ++initval;
    assert(PointBuff_size(pb) == 0);

    PointBuff_push(pb, cx, npar, x);
    assert(PointBuff_size(pb) == 1);

    PointBuff_push(pb, cy, npar, y);
    assert(PointBuff_size(pb) == 2);

    PointBuff_push(pb, cz, npar, z);
    assert(PointBuff_size(pb) == 2);

    double c, a[npar];

    c = PointBuff_pop(pb, npar, a);
    assert(PointBuff_size(pb) == 1);
    assert(c == cy);
    assert(0 == memcmp(a, y, sizeof(a)));

    c = PointBuff_pop(pb, npar, a);
    assert(PointBuff_size(pb) == 0);
    assert(c == cz);
    assert(0 == memcmp(a, z, sizeof(a)));

    PointBuff_free(pb);
    unitTstResult("PointBuff", "OK");
}
