#ifndef LEGO_COST
#define LEGO_COST

struct ThreadData {
    unsigned    seed;
};

/// Parameters of cost function--that which is minimized.
typedef struct CostPar {
    const BranchTab *obs;      // observed site pattern frequencies; normalized
    const GPTree *gptree;   // model of population history
    int         nThreads; // number of threads to use
    long        nreps;    // number of repetitions to simulate
} CostPar;

double      costFun(int dim, double x[dim], void *jdata, void *tdata);

#endif
