#ifndef BRANCHTAB_INCLUDED
#  define BRANCHTAB_INCLUDED

#include "typedefs.h"

BranchTab    *BranchTab_new(void);
void          BranchTab_free(BranchTab * self);
double        BranchTab_get(BranchTab * self, tipId_t tipid);
double        BranchTab_tabulate(BranchTab * self, tipId_t tipid,
                                 double branch);
unsigned long BranchTab_size(BranchTab * self);
void          BranchTab_print(BranchTab *self);
#endif
