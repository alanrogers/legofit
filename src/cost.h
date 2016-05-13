#ifndef LEGO_COST
#define LEGO_COST

/// Parameters of cost function--that which is minimized.
typedef struct CostPar {
    int         npat
    double     *obs;
    GPTree     *gptree;
} CostPar;

double      costFun(int dim, double x[dim], void *jdata, void *tdata);

#endif
