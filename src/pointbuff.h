#ifndef ARR_POINTBUFF_H
#define ARR_POINTBUFF_H

PointBuff *PointBuff_new(unsigned npar, unsigned totpts);
void     PointBuff_free(Wraparound * self);
unsigned PointBuff_size(const PointBuff * self);
void     PointBuff_push(PointBuff * self, double cost, int n, double param[n]);
double   PointBuff_pop(PointBuff * self, int n, double param[n]);

#endif
