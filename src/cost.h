#ifndef LEGO_COST
#define LEGO_COST

#include "typedefs.h"

struct ThreadData {
    unsigned    seed;
};

/// Parameters of cost function--that which is minimized.
typedef struct CostPar {
    BranchTab *obs;       // observed site pattern frequencies; normalized
    GPTree     *gptree;   // model of population history
    int         nThreads; // number of threads to use
    int         doSing;   // nonzero => use singleton site patterns
#if COST!=KL_COST
    double      u;        // mutation rate per generation
    long        nnuc;     // number of nucleotide sites in genome
#endif
    SimSched   *simSched;
} CostPar;

double      costFun(int dim, double x[dim], void *jdata, void *tdata);
void       *CostPar_dup(const void * arg);
void        CostPar_free(void *arg);

#endif
