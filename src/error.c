#include "error.h"
#include <stdio.h>
#include <string.h>

// Returns 0 on success, 1 on buffer overflow.
int mystrerror_r(int errnum, char *buff, size_t len) {
    int status=0, rval=0;

    switch(errnum) {
    case NO_ANCESTRAL_ALLELE:
        rval = snprintf(buff, len, "No ancestral allele");
        break;
    case MONOMORPHIC_SITE:
        rval = snprintf(buff, len, "Monomorphic site");
        break;
    case REF_MISMATCH:
        rval = snprintf(buff, len, "Inconsistent REF alleles");
        break;
    case ALLELE_MISMATCH:
        rval = snprintf(buff, len, "Inconsistent alleles");
        break;
    case MULTIPLE_ALT:
        rval = snprintf(buff, len, "Multiple ALT alleles");
        break;
    case BUFFER_OVERFLOW:
        rval = snprintf(buff, len, "Buffer overflow");
        break;
    case BAD_RAF_INPUT:
        rval = snprintf(buff, len, "Bad .raf input file");
        break;
    case BAD_DAF_INPUT:
        rval = snprintf(buff, len, "Bad .daf input file");
        break;
    case BAD_SORT:
        rval = snprintf(buff, len, "Incorrect sort");
        break;
    case DIMENSION_MISMATCH:
        rval = snprintf(buff, len, "Inconsistent array dimensions");
        break;
    case NAME_MISMATCH:
        rval = snprintf(buff, len, "Inconsistent variable names");
        break;
    case TOO_MANY_PARENTS:
        rval = snprintf(buff, len, "Too many parents");
        break;
    case TOO_MANY_CHILDREN:
        rval = snprintf(buff, len, "Too many children");
        break;
    case DATE_MISMATCH:
        rval = snprintf(buff, len, "Date mismatch");
        break;
    default:
        status = strerror_r(errnum, buff, len);
    }
    if(rval>=len || status!=0) {
        buff[len-1] = '\0';
        fprintf(stderr,"%s:%s:%d: buffer overflow\n",
                __FILE__,__func__,__LINE__);
        return 1;
    }

    return 0;
}
