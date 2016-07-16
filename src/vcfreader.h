#ifndef VCFREADER_INCLUDED
#define VCFREADER_INCLUDED

#include "typedefs.h"

VCFReader *VCFReader_new(const char *fname);
void VCFReader_free(VCFReader *self);
void VCFReader_parseHdr(VCFReader *self);
int VCFReader_next(VCFReader *self);

#endif
