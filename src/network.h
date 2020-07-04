#ifndef ARR_NETWORK_H
#  define ARR_NETWORK_H

#  include "typedefs.h"
#  include <stdio.h>
#  include <gsl/gsl_rng.h>

#ifndef EXTERN
#  define EXTERN extern
#endif

EXTERN enum NetworkType networkType;

EXTERN void       *(*Network_dup)(const void *old);
EXTERN int         (*Network_feasible)(const void *self, int verbose);
EXTERN void        (*Network_free)(void *self);
EXTERN LblNdx      (*Network_getLblNdx)(void *self);
EXTERN const char *(*Network_getNameFree)(void * self, int i);
EXTERN void        (*Network_getParams)(void *self, int n, double x[n]);
EXTERN void       *(*Network_new)(const char *fname, Bounds bnd);
EXTERN int         (*Network_nFree)(const void *self);
EXTERN void        (*Network_patprob)(void *self, BranchTab *branchtab,
                                      gsl_rng *rng, unsigned long nreps,
                                      int doSing);
EXTERN void        (*Network_printParStore)(void *self, FILE *fp);
EXTERN void        (*Network_printParStoreFree)(void *self, FILE *fp);
EXTERN void        (*Network_randomize)(void *self, gsl_rng *rng);
EXTERN void        (*Network_sanityCheck)(void *self, const char *file,
                                          int line);
EXTERN int         (*Network_setParams)(void *self, int n, double x[n]);
EXTERN void        (*Network_initStateVec)(void *gpt, int ndx, int n,
                                           double x[n], gsl_rng *rng);

EXTERN void       *(*Node_new)(int twoN_i, int start_i, ParStore *ps,
                               NodeStore *ns);
EXTERN int         (*Node_addChild)(void * parent, void * child);
EXTERN int         (*Node_mix)(void * child, int mix_i, void * introgressor,
                               void * native, ParStore *ps);
EXTERN void       *(*Node_root)(void * self);
EXTERN void        (*Node_print)(FILE * fp, void * vself, int indent);

void Network_init(enum NetworkType type);

#endif


