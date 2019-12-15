#ifndef ARR_STIRLING2_H
#define ARR_STIRLING2_H

typedef struct Stirling2 Stirling2;

void Stirling2_free(Stirling2 *self);
Stirling2 *Stirling2_new(long unsigned nmax);
long unsigned Stirling2_val(Stirling2 *self, long unsigned n, long unsigned k);

#endif
