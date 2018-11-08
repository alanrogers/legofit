#ifndef ARR_GETLINETOK_H
#define ARR_GETLINETOK_H

#include "typedefs.h"
#include <stdio.h>

GetLineTok *GetLineTok_new(size_t buffsize, int maxtokens, FILE *fp);
int GetLineTok_next(GetLineTok *self, const char *sep,
                        const char *extraneous);
int GetLineTok_ntokens(GetLineTok *self);
char *GetLineTok_token(GetLineTok * self, int index);
void GetLineTok_print(GetLineTok *self, FILE *ofp);
void GetLineTok_free(GetLineTok *self);

#endif
