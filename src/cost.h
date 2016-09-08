#ifndef LEGO_COST
#define LEGO_COST

#include "typedefs.h"

struct ThreadData {
    unsigned    seed;
};

/// Parameters of cost function--that which is minimized.
typedef struct CostPar {
    const BranchTab *obs;      // observed site pattern frequencies; normalized
    GPTree     *gptree;   // model of population history
    int         nThreads; // number of threads to use
    long        nreps;    // number of repetitions to simulate
    int         doSing;   // nonzero => use singleton site patterns
    double      u;        // mutation rate per generation
    long        nnuc;     // number of nucleotide sites in genome
} CostPar;

double      costFun(int dim, double x[dim], void *jdata, void *tdata);

#endif
