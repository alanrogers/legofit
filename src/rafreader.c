/**
   @file rafreader.c
   @brief Class RAFReader: read a raf file.

   @copyright Copyright (c) 2016, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "rafreader.h"
#include "tokenizer.h"
#include "misc.h"
#include "error.h"
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <unistd.h>

#define MAXFIELDS 5

int iscomment(const char *s);
int RAFReader_cmp(const RAFReader *lhs, const RAFReader *rhs);
void RAFReader_popen(RAFReader *self);

/// RAFReader constructor
RAFReader *RAFReader_new(const char *fname) {
    RAFReader *self = malloc(sizeof(*self));
    CHECKMEM(self);
    memset(self, 0, sizeof(RAFReader));
    self->fname = strdup(fname);
    CHECKMEM(self->fname);

    // Does input file exist?
    if( access( self->fname, F_OK ) == -1 ) {
        fprintf(stderr, "%s:%d: nonexistent input:  %s\n",
                __FILE__,__LINE__, self->fname);
        exit(EXIT_FAILURE);
    }

    // Is input file compressed?
    const char *gz = ".gz";
    char *pos = strstr(self->fname, gz);
    if(pos == gz) {
        fprintf(stderr,"%s:%d: input file name \"%s\" is empty\n",
                __FILE__,__LINE__, self->fname);
        exit(EXIT_FAILURE);
    }else if(pos == NULL) {  // not compressed
        self->fp = fopen(self->fname, "r");
        self->ispipe = 0;
    }else {                  // compressed
        RAFReader_popen(self);
    }
    if(self->fp == NULL) {
        fprintf(stderr, "%s:%s:%d: can't open \"%s\" for input.\n",
                __FILE__, __func__, __LINE__, self->fname);
        exit(EXIT_FAILURE);
    }

    self->tkz = Tokenizer_new(MAXFIELDS);
    self->snpid = -1;
    self->raf = strtod("NaN", NULL);
    self->daf = strtod("NaN", NULL);
    return self;
}

// open pipe to read compressed file
void RAFReader_popen(RAFReader *self) {
    
    self->ispipe = 1;
    char cmd[10000];
    int status = snprintf(cmd, sizeof cmd, "gunzip -c %s", self->fname);
    if(status >= sizeof cmd) {
        fprintf(stderr, "%s:%d: ERR: Input filename is too large."
                " Max: %zu\n",
                __FILE__, __LINE__, sizeof(cmd) - 1);
        exit(EXIT_FAILURE);
    }
#ifdef _WIN32
    // windows
    self->fp = _popen(cmd, "r");
#else
    // osx, linux, or unix
    self->fp = popen(cmd, "r");
#endif
    if(self->fp == NULL) {
        fprintf(stderr,"%s:%s:%d: can't open pipe with cmd %s\n",
                __FILE__,__func__,__LINE__, cmd);
        exit(EXIT_FAILURE);
    }
}

/// Clear all chromosome names
void RAFReader_clearChromosomes(int n, RAFReader * r[n]) {
    int i;
    for(i = 0; i < n; ++i)
        r[i]->chr[0] = '\0';
}

/// RAFReader destructor
void RAFReader_free(RAFReader * self) {
    if(self->ispipe) {
#ifdef _WIN32
        _pclose(self->fp);
#else
        pclose(self->fp);
#endif
    }else
        fclose(self->fp);
    free(self->fname);
    Tokenizer_free(self->tkz);
    free(self);
}

/// Return 1 if first non-white character in string is '#'; 0
/// otherwise.
int iscomment(const char *s) {
    int rval;
    while(*s != '\0' && isspace(*s))
        ++s;
    rval = (*s == '#');
    return rval;
}

/// Read the next site and set derived allele frequency (daf) within
/// each reader.  @return 0 on success, or else EOF,
/// NO_ANCESTRAL_ALLELE, BUFFER_OVERFLOW, BAD_RAF_INPUT, BAD_SORT,
/// MONOMORPHIC_SITE, or an errno code for failure to parse a
/// floating-point number.
int RAFReader_next(RAFReader * self) {
    int ntokens;
    int status;
    char buff[1024];
    long unsigned prevnucpos = 0UL;

    // Find a line of input
    while(1) {
        if(fgets(buff, sizeof(buff), self->fp) == NULL)
            return EOF;
        if(NULL == strchr(buff, '\n') && !feof(self->fp)) {
#ifdef NDEBUG
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
#endif
            return BUFFER_OVERFLOW;
        }
        if(iscomment(buff))
            continue;
        Tokenizer_split(self->tkz, buff, "\t");
        ntokens = Tokenizer_strip(self->tkz, " \n");
        if( ntokens == 5) {
            // reference allele must be a single nucleotide
            if(1 != strlen(Tokenizer_token(self->tkz, 2)))
                continue;
            // alternate allele must be a single nucleotide
            if(1 != strlen(Tokenizer_token(self->tkz, 3)))
                continue;
        }
        if(ntokens > 0)
            break;
    }

    if(ntokens != 5) {
#ifdef NDEBUG
        fprintf(stderr, "%s:%d: Each line of .raf file must have 5 tokens,"
                " but current line has %d.\n", __FILE__, __LINE__, ntokens);
#endif
        return BAD_RAF_INPUT;
    }

    ++self->snpid;

    // Chromosome
    char prev[RAFSTRSIZE];
    assert(sizeof prev == sizeof self->chr);
    memcpy(prev, self->chr, sizeof prev);
    status = snprintf(self->chr, sizeof self->chr, "%s",
                      Tokenizer_token(self->tkz, 0));
    if(status >= sizeof self->chr) {
#ifdef NDEBUG
        fprintf(stderr, "%s:%d: chromosome name too long: %s\n",
                __FILE__, __LINE__, Tokenizer_token(self->tkz, 0));
#endif
        return BUFFER_OVERFLOW;
    }
    int diff = strcmp(prev, self->chr);
    if(diff > 0) {
#ifdef NDEBUG
        fprintf(stderr, "%s:%s:%d: Chromosomes missorted in input.\n",
                __FILE__, __func__, __LINE__);
        fprintf(stderr, "          \"%s\" precedes \"%s\".\n",
                prev, self->chr);
        Tokenizer_print(self->tkz, stderr);
#endif
        return BAD_SORT;
    } else if(diff < 0) {
        // new chromosome
        prevnucpos = 0UL;
    } else
        prevnucpos = self->nucpos;

    // Nucleotide position
    self->nucpos = strtoul(Tokenizer_token(self->tkz, 1), NULL, 10);
    if(prevnucpos == self->nucpos) {
#ifdef NDEBUG
        fprintf(stderr, "%s:%d: Duplicate line in raf file. chr=%s pos=%lu\n",
                __FILE__, __LINE__, self->chr, self->nucpos);
#endif
        return BAD_SORT;
    } else if(prevnucpos > self->nucpos) {
#ifdef NDEBUG
        fprintf(stderr, "%s:%d: positions missorted chr=%s "
                "prev=%lu curr=%lu\n",
                __FILE__, __LINE__, self->chr, prevnucpos, self->nucpos);
#endif
        return BAD_SORT;
    }
    // Reference allele
    status = snprintf(self->ref, sizeof(self->ref), "%s",
                      Tokenizer_token(self->tkz, 2));
    strlowercase(self->ref);

    // Alternate allele
    snprintf(self->alt, sizeof(self->alt), "%s",
             Tokenizer_token(self->tkz, 3));
    strlowercase(self->alt);

    // Reference allele frequency
    char *token, *end;
    token = Tokenizer_token(self->tkz, 4);
    errno = 0;
    self->raf = strtod(token, &end);
    if(end == token)
        errno = EINVAL;
    if(errno) {
        char err_buff[50];
        strerror_r(errno, err_buff, sizeof(err_buff));
#ifdef NDEBUG
        fprintf(stderr, "%s:%d: Bad float \"%s\" (%s); chr=%s pos=%lu\n",
                __FILE__, __LINE__, token, err_buff, self->chr, self->nucpos);
#endif
        return errno;
    }

    // Derived allele frequency is set later
    self->daf = strtod("NaN", NULL);

    return 0;
}

/// Rewind raf file.
/// @return 0 on success; -1 on failure
int RAFReader_rewind(RAFReader * self) {
    if(self->ispipe) {
#ifdef _WIN32
        _pclose(self->fp);
#else
        pclose(self->fp);
#endif
        RAFReader_popen(self);
        if(self->fp == NULL)
            return -1;
        return 0;
    }
    return fseek(self->fp, 0L, SEEK_SET);
}

/**
   Compare two RAFReader objects. If the "chr" fields differ, return
   positive if lhs->chr > rhs->chr; return negative if the reverse
   inequality holds. Otherwise return positive, 0, or negative to
   match the sign of lhs->nucpos - rhs->nucpos.
 */
int RAFReader_cmp(const RAFReader *lhs, const RAFReader *rhs) {
    int c = strcmp(lhs->chr, rhs->chr);
    if(c)
        return c;
    if(lhs->nucpos > rhs->nucpos)
        return 1;
    if(lhs->nucpos < rhs->nucpos)
        return -1;
    return 0;
}

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) > (Y) ? (Y) : (X))

/// Advance an array of RAFReaders to the next shared position, and
/// set derived allele frequency within each RAFReader.
/// @param[in] n number of RAFReader objects in array
/// @param[in] r array of RAFReader objects. Last one should be outgroup.
/// @return 0 on success, or one of several error codes on failure.
int RAFReader_multiNext(int n, RAFReader * r[n]) {
    int i, imax, status, diff, atsamepos;

    assert(n <= 8*sizeof(bits_t));

    for(i = 0; i < n; ++i) {
        if( (status = RAFReader_next(r[i])) )
            return status;
    }

    // Bits indicate whether readers are at maximum position.
    // We begin by assuming that reader 0 is the maximum, so
    // the initial value of atmax has only the 0th bit turned
    // on.
    const bits_t unity = 1u;
    bits_t atmax = unity;
    imax = 0;

    while(1) {

        bits_t currbit;

        // Set values of atmax and atsamepos. The i'th bit of
        // atmax is 1 if the i'th reader is at the maximum position
        // and is 0 otherwise. atsamepos is 1 if all readers are
        // at the same position, 0 otherwise.
        atsamepos = 1;
        for(i = 0; i < n; ++i) {
            currbit = unity << i;  // bit i on; all others off
            if( atmax & currbit ) {
                // reader i is at current maximum
                continue;
            }
            diff = RAFReader_cmp(r[i], r[imax]);
            if(diff > 0) {
                // found new maximum
                atmax = currbit;
                imax = i;
                atsamepos = 0;
            } else if(diff < 0) {
                atsamepos = 0;
            }else
                atmax |= currbit;  // turn on bit i in atmax
        }

        // If all readers are at the same position, then we're done.
        if(atsamepos) {
            assert(n == num1bits(atmax));
            break;
        }

        // Increment readers that are not at the maximum position.
        for(i=0; i < n; ++i) {
            currbit = unity << i;

            while(0 == (atmax & currbit) ) {
                if( (status = RAFReader_next(r[i])) )
                    return status;
                diff = RAFReader_cmp(r[i], r[imax]);
                if(diff > 0) {
                    // found new maximum
                    atmax = currbit;
                    imax = i;
                } else if(diff == 0) {
                    atmax |= currbit;
                }
            }
        }
    }

    // Make sure REF and ALT are consistent across readers
    if((status = RAFReader_alleleCheck(n, r)))
        return status;

    return 0;
}

/// Set derived allele frequency within each RAFReader.
/// @param[in] n number of RAFReader objects in array
/// @param[in] r array of RAFReader objects. Last one should be outgroup.
/// @return 0 on success, or one of several error codes on failure.
int RAFReader_findDaf(int n, RAFReader * r[n]) {

    int i, status;
    
    // Make sure reference allele isn't fixed in readers, excluding
    // the outgroup (reader n-1). If it's fixed, then we can't call
    // the ancestral allele. 
    double maxp, minp;
    maxp = minp = RAFReader_raf(r[0]);
    for(i = 1; i < n - 1; ++i) {
        double p = RAFReader_raf(r[i]); // reference allele freq
        minp = fmin(minp, p);
        maxp = fmax(maxp, p);
    }
    if(maxp == 0.0 || minp == 1.0)
        return MONOMORPHIC_SITE;

    // Make sure REF and ALT are consistent across readers
    if((status = RAFReader_alleleCheck(n, r)))
        return status;

    // Set derived allele frequency
    double ogf = r[n - 1]->raf; // freq of ref in outgroup
    if(ogf == 0.0) {
        // reference allele is derived
        for(i = 0; i < n; ++i)
            r[i]->daf = r[i]->raf;
    } else if(ogf == 1.0) {
        // alternate allele is derived
        for(i = 0; i < n; ++i)
            r[i]->daf = 1.0 - r[i]->raf;
    } else {
        // outgroup is polymorphic: can't call ancestral allele
        return NO_ANCESTRAL_ALLELE;
    }

    return 0;
}


/// Return 0 if ref and alt alleles of all readers match; return
/// REF_MISMATCH if there is a mismatch in REF alleles; return
/// MULTIPLE_ALT if there is a mismatch in ALT alleles.
int RAFReader_alleleCheck(int n, RAFReader * r[n]) {
    char *ref = r[0]->ref;
    char *alt = r[0]->alt;
    int altMissing = (0 == strcmp(".", alt));
    int i;
    for(i = 1; i < n; ++i) {
        if(0 != strcmp(ref, r[i]->ref))
            return REF_MISMATCH;
        int currAltMissing = (0 == strcmp(".", r[i]->alt));
        if(altMissing && !currAltMissing) {
            altMissing = 0;
            alt = r[i]->alt;
            strcpy(r[0]->alt, alt);
            continue;
        }
        if(!altMissing && !currAltMissing && 0 != strcmp(alt, r[i]->alt))
            return MULTIPLE_ALT;
    }
    return 0;
}

/// Print header for raf file.
void RAFReader_printHdr(FILE * fp) {
    fprintf(fp, "%30s %5s %10s %3s %3s %8s %8s\n",
            "file", "chr", "pos", "ref", "alt", "raf", "daf");
}

/// Print current line of raf file
void RAFReader_print(RAFReader * r, FILE * fp) {
    assert(r->fname);
    int status;
    char buff[500];
    status = snprintf(buff, sizeof buff, "%s", r->fname);
    if(status >= sizeof buff) {
        fprintf(stderr, "%s:%s:%d: buffer overflow\n",
                __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "%30s %5s %10lu %3s %3s %8.6lg %8.6lg\n",
            strltrunc(buff, 30), r->chr, r->nucpos, r->ref, r->alt, r->raf,
            r->daf);
}

void RAFReader_printArray(int n, RAFReader * r[n], FILE * fp) {
    int i;
    for(i = 0; i < n; ++i)
        RAFReader_print(r[i], fp);
}

/// Return reference allele frequency of current line of raf file.
double RAFReader_raf(RAFReader * r) {
    return r->raf;
}

/// Return derived allele frequency of current line of raf file.
double RAFReader_daf(RAFReader * r) {
    return r->daf;
}
