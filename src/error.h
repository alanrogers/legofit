#ifndef ERROR_INCLUDED
#define ERROR_INCLUDED

#include <stdlib.h>

// Number from 1000 to avoid conflict with errno
enum {NO_ANCESTRAL_ALLELE=1000,
      REF_MISMATCH,
      ALLELE_MISMATCH,
      MULTIPLE_ALT,
      BUFFER_OVERFLOW,
      BAD_RAF_INPUT,
      BAD_DAF_INPUT,
      BAD_SORT};

int mystrerror_r(int errnum, char *buff, size_t len);

#endif
