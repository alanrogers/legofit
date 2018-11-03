#ifndef ARR_GETLINETOK_H
#define ARR_GETLINETOK_H

#include "typedefs.h"

GetLineTok *GetLineTok_new(size_t buffsize, int maxtokens, FILE *fp);
int GetLineTok_nextLine(GetLineTok *self, const char *sep,
                        const char *extraneous);
int GetLineTok_ntokens(GetLineTok *self);
char *GetLineTok_token(Tokenizer * t, int index);


#endif
