#ifdef ERROR_INCLUDED
#define ERROR_INCLUDED

// Number from 1000 to avoid conflict with errno
enum {NO_ANCESTRAL_ALLELE=1000,
      REF_ALT_MISMATCH,
      BUFFER_OVERFLOW,
      BAD_RAF_INPUT,
      BAD_SORT};

int mystrerror_r(int errnum, char *buff, size_t len);

#endif
