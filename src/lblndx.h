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
const char *LblNdx_lbl(const LblNdx * self, unsigned i);
unsigned    LblNdx_size(const LblNdx * self);
void        LblNdx_sanityCheck(const LblNdx *self, const char *file, int line);
int         LblNdx_equals(const LblNdx *lhs, const LblNdx *rhs);
tipId_t     LblNdx_getTipId(const LblNdx *self, const char *lbl);
void        LblNdx_print(const LblNdx *self, FILE *fp);

char       *patLbl(size_t n, char buff[n], tipId_t tid, const LblNdx * lblndx);

#endif
