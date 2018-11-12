#ifndef ARR_LINEREADER_H
#define ARR_LINEREADER_H

#include "typedefs.h"
#include <stdio.h>

LineReader *LineReader_new(size_t buffsize);
char *LineReader_next(LineReader *self, FILE *fp);
void LineReader_free(LineReader *self);

#endif
