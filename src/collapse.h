#ifndef ARR_COLLAPSE_H
#define ARR_COLLAPSE_H

#include "typedefs.h"
#include "lblndx.h"

struct LblNdxBranchTab {
    LblNdx lndx;
    BranchTab *branchtab;
};

LblNdxBranchTab removePops(LblNdx lndx, BranchTab *bt, const char *deleteStr);
LblNdxBranchTab collapsePops(LblNdx lndx, BranchTab *bt,
                             const char *collapseStr,
                             const char *lbl);
#endif
