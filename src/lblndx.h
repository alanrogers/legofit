#ifndef ARR_LBLNDX
#  define ARR_LBLNDX

#  include "gptree.h"
#  include "popnode.h"

struct LblNdx {
    unsigned    n;              // number of samples
    char        lbl[MAXSAMP][POPNAMESIZE];
};

void        LblNdx_init(LblNdx * self);
void        LblNdx_addSamples(LblNdx * self, unsigned nsamples,
							   const char *lbl);
const char *LblNdx_lbl(LblNdx * self, unsigned i);
unsigned    LblNdx_size(LblNdx * self);
void        LblNdx_sanityCheck(LblNdx *self, const char *file, int line);

int         LblNdx_equals(LblNdx *lhs, LblNdx *rhs);

#endif
