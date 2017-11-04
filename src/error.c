int mystrerror_r(int errnum, char *buff, size_t len) {
    int status=0, rval;
    
    switch(errnum) {
    case NO_ANCESTRAL_ALLELE:
        rval = snprintf(buff, len, "No ancestral allele");
        if(rval >= buflen)
            status BUFFER_OVERFLOW;
        break;
    case REF_ALT_MISMATCH:
        rval = snprintf(buff, len, "Inconsistent REF and ALT alleles");
        if(rval >= buflen)
            status BUFFER_OVERFLOW;
        break;
    case BUFFER_OVERFLOW:
        rval = snprintf(buff, len, "Buffer overflow");
        if(rval >= buflen)
            status BUFFER_OVERFLOW;
        break;
    case BAD_RAF_INPUT:
        rval = snprintf(buff, len, "Bad .raf input file");
        if(rval >= buflen)
            status BUFFER_OVERFLOW;
        break;
    case BAD_SORT:
        rval = snprintf(buff, len, "Incorrect sort");
        if(rval >= buflen)
            status BUFFER_OVERFLOW;
        break;
    default:
        status = strerror_r(errnum, buff, len);
    }
    return status;
}
