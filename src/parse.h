#ifndef HAVE_PARSE
#define HAVE_PARSE

#include "typedefs.h"

PopNode    *parse(FILE * fp, HashTab * ht);
PopNode    *mktree(FILE * fp, HashTab * ht);

#endif
