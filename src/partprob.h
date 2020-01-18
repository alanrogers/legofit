#ifndef PARTPROB_H
#define PARTPROB_H

long double lnCoalConst(unsigned n, unsigned k);
double probPartition(unsigned k, unsigned y[k], long double lnconst);

#endif
