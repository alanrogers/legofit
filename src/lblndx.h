#ifndef ARR_LBLNDX
#  define ARR_LBLNDX

#  include "gptree.h"
#  include "popnode.h"
#  include "strdblqueue.h"

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
int         LblNdx_collapse(LblNdx *self, tipId_t collapse, const char *lbl);
int         LblNdx_rmPops(LblNdx *self, tipId_t remove);
int         LblNdx_from_StrDblQueue(LblNdx *lndx, StrDblQueue *queue);
char       *patLbl(size_t n, char buff[n], tipId_t tid, const LblNdx * lblndx);
void        orderpat(int n, unsigned order[n], tipId_t tid[n]);
int         compare_tipId(const void *void_x, const void *void_y);

#endif
