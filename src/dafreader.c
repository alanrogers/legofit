/**
   @file dafreader.c
   @brief Class DAFReader: read a daf file.

   @copyright Copyright (c) 2016, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "dafreader.h"
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
int DAFReader_cmp(const DAFReader *lhs, const DAFReader *rhs);

/// DAFReader constructor
DAFReader  *DAFReader_new(const char *fname) {
    DAFReader  *self = malloc(sizeof(*self));
    CHECKMEM(self);
    memset(self, 0, sizeof(DAFReader));
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
    }else if(pos == NULL) {
        // not compressed
        self->fp = fopen(self->fname, "r");
        self->ispipe = 0;
    }else {
        // compressed
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
    }
    if(self->fp == NULL) {
        fprintf(stderr, "%s:%s:%d: can't open \"%s\" for input.\n",
                __FILE__, __func__, __LINE__, self->fname);
        exit(EXIT_FAILURE);
    }

    self->tkz = Tokenizer_new(MAXFIELDS);
    self->snpid = -1;
    self->p = strtod("NaN", NULL);
    return self;
}

/// Clear all chromosome names
void DAFReader_clearChromosomes(int n, DAFReader * r[n]) {
    int         i;
    for(i = 0; i < n; ++i)
        r[i]->chr[0] = '\0';
}

/// DAFReader destructor
void DAFReader_free(DAFReader * self) {
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
    int         rval;
    while(*s != '\0' && isspace(*s))
        ++s;
    rval = (*s == '#');
    return rval;
}

/// Read the next site.
/// @return 0 on success; EOF on end of file; abort on other errors.
int DAFReader_next(DAFReader * self) {
    int         ntokens1;
    int         ntokens;
    int         status;
    char        buff[1024];
    long unsigned prevnucpos = 0UL;

    // Find a line of input
    while(1) {
        if(fgets(buff, sizeof(buff), self->fp) == NULL)
            return EOF;
        if(NULL == strchr(buff, '\n') && !feof(self->fp)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
            return BUFFER_OVERFLOW;
        }
        if(iscomment(buff))
            continue;
        ntokens1 = Tokenizer_split(self->tkz, buff, " ");
        ntokens = Tokenizer_strip(self->tkz, " \n");
        if( ntokens == 5) {
            // ancestral allele must be a single nucleotide
            if(1 != strlen(Tokenizer_token(self->tkz, 2)))
                continue;
            // derived allele must be a single nucleotide
            if(1 != strlen(Tokenizer_token(self->tkz, 3)))
                continue;
        }
        if(ntokens > 0)
            break;
    }

    if(ntokens != 5) {
        fprintf(stderr, "%s:%d: Each line of .daf file must have 5 tokens,"
                " but current line has %d.\n", __FILE__, __LINE__, ntokens);
        fprintf(stderr, "ntokens1=%d\n", ntokens1);
        fprintf(stderr, "buff: %s\n", buff);
        Tokenizer_printSummary(self->tkz, stderr);
        Tokenizer_print(self->tkz, stderr);
        return BAD_DAF_INPUT;
    }

    ++self->snpid;

    // Chromosome
    char        prev[DAFSTRSIZE];
    assert(sizeof prev == sizeof self->chr);
    memcpy(prev, self->chr, sizeof prev);
    status = snprintf(self->chr, sizeof self->chr, "%s",
                      Tokenizer_token(self->tkz, 0));
    if(status >= sizeof self->chr) {
        fprintf(stderr, "%s:%d: chromosome name too long: %s\n",
                __FILE__, __LINE__, Tokenizer_token(self->tkz, 0));
        return BUFFER_OVERFLOW;
    }
    int         diff = strcmp(prev, self->chr);
    if(diff > 0) {
        fprintf(stderr, "%s:%s:%d: Chromosomes missorted in input.\n",
                __FILE__, __func__, __LINE__);
        fprintf(stderr, "          \"%s\" precedes \"%s\".\n",
                prev, self->chr);
        Tokenizer_printSummary(self->tkz, stderr);
        Tokenizer_print(self->tkz, stderr);
        return BAD_SORT;
    } else if(diff < 0) {
        // new chromosome
        prevnucpos = 0UL;
    } else
        prevnucpos = self->nucpos;

    // Nucleotide position
    self->nucpos = strtoul(Tokenizer_token(self->tkz, 1), NULL, 10);
    if(prevnucpos == self->nucpos) {
        fprintf(stderr, "%s:%d: Duplicate line in daf file. chr=%s pos=%lu\n",
                __FILE__, __LINE__, self->chr, self->nucpos);
        return BAD_SORT;
    } else if(prevnucpos > self->nucpos) {
        fprintf(stderr, "%s:%d: positions missorted chr=%s "
                "prev=%lu curr=%lu\n",
                __FILE__, __LINE__, self->chr, prevnucpos, self->nucpos);
        return BAD_SORT;
    }
    // Ancestral allele
    status = snprintf(self->aa, sizeof(self->aa), "%s",
                      Tokenizer_token(self->tkz, 2));
    strlowercase(self->aa);

    // Derived allele
    snprintf(self->da, sizeof(self->da), "%s", Tokenizer_token(self->tkz, 3));
    strlowercase(self->da);

    // Derived allele frequency
    char *token, *end;
    token = Tokenizer_token(self->tkz, 4);
    errno=0;
    self->p = strtod(token, &end);
    if(end==token)
        errno = EINVAL;
    if(errno) {
        char err_buff[500];
        strerror_r(errno, err_buff, sizeof(err_buff));
        fprintf(stderr,"%s:%d: Bad float \"%s\" (%s); chr=%s pos=%lu\n",
                __FILE__,__LINE__, token, err_buff, self->chr, self->nucpos);
        return errno;
    }

    return 0;
}

/// Rewind daf file.
/// @return 0 on success; -1 on failure
int DAFReader_rewind(DAFReader * self) {
    return fseek(self->fp, 0L, SEEK_SET);
}

/**
   Compare two DAFReader objects. If the "chr" fields differ, return
   positive if lhs->chr > rhs->chr; return negative if the reverse
   inequality holds. Otherwise return positive, 0, or negative to
   match the sign of lhs->nucpos - rhs->nucpos.
 */
int DAFReader_cmp(const DAFReader *lhs, const DAFReader *rhs) {
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

#if 1

int DAFReader_multiNext(int n, DAFReader * r[n]) {
    int i, imax, status, diff, atsamepos;

    for(i = 0; i < n; ++i) {
        if( (status = DAFReader_next(r[i])) )
            return status;
    }

    // bits indicate whether readers are at maximum
    const tipId_t unity = 1u;
    tipId_t atmax = unity;
    imax = 0;

    while(1) {

        tipId_t currbit;

        // Set values of atmax and atsamepos. The i'th bit of
        // atmax is 1 if the i'th reader is at the maximum position
        // and is 0 otherwise. atsamepos is 1 if all readers are
        // at the same position, 0 otherwise.
        atsamepos = 1;
        for(i = 0; i < n; ++i) {
            currbit = unity << i;
            if( atmax & currbit ) {
                // reader i is at current maximum
                continue;
            }
            diff = DAFReader_cmp(r[i], r[imax]);
            if(diff > 0) {
                // found new maximum
                atmax = currbit;
                imax = i;
                atsamepos = 0;
            } else if(diff < 0) {
                atsamepos = 0;
            }else
                atmax |= currbit;
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
                if( (status = DAFReader_next(r[i])) )
                    return status;
                diff = DAFReader_cmp(r[i], r[imax]);
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

    // Make sure reference allele isn't fixed in readers. If it's
    // fixed, then we can't call the ancestral allele.
    double maxp, minp;
    maxp = minp = DAFReader_daf(r[0]);
    for(i=1; i < n; ++i) {
        double p = DAFReader_daf(r[i]); // derived allele freq
        minp = fmin(minp, p);
        maxp = fmax(maxp, p);
    }
    if(maxp == 0.0 || minp == 1.0)
        return MONOMORPHIC_SITE;

    if(!DAFReader_allelesMatch(n, r))
        return ALLELE_MISMATCH;

    return 0;
}

#else

/// Advance an array of DAFReaders to the next shared position.
/// @return 0 on success or EOF on end of file.
int DAFReader_multiNext(int n, DAFReader * r[n]) {
    int   i, status;
    unsigned long maxnuc = 0, minnuc = ULONG_MAX;
    int   imaxchr;   // index of reader with maximum chromosome position
    int   onSameChr; // indicates whether all readers are on same chromosome.
    int   diff;
    char  currchr[DAFSTRSIZE] = { '\0' }; // current chromosome

    // Set index, imaxchr, of reader with maximum
    // chromosome values in lexical sort order, and
    // set boolean flag, onSameChr, which indicates
    // whether all readers are on same chromosome.
    if( (status = DAFReader_next(r[0])) )
        return status;

    imaxchr = 0;
    onSameChr = 1;
    for(i = 1; i < n; ++i) {
        if( (status = DAFReader_next(r[i])) )
            return status;

        diff = strcmp(r[i]->chr, r[imaxchr]->chr);
        if(diff > 0) {
            onSameChr = 0;
            imaxchr = i;
        } else if(diff < 0)
            onSameChr = 0;
    }

    // Loop until both chr and position are homogeneous.
    do {
        // get them all on the same chromosome
        while(!onSameChr) {
            onSameChr = 1;
            for(i = 0; i < n; ++i) {
                if(i == imaxchr)
                    continue;
                while((diff = strcmp(r[i]->chr, r[imaxchr]->chr)) < 0) {
                    if( (status = DAFReader_next(r[i])) )
                        return status;
                }
                assert(diff >= 0);
                if(diff > 0) {
                    imaxchr = i;
                    onSameChr = 0;
                }
            }
        }

        assert(onSameChr);
        maxnuc = minnuc = r[0]->nucpos;
        for(i = 1; i < n; ++i) {
            maxnuc = MAX(maxnuc, r[i]->nucpos);
            minnuc = MIN(minnuc, r[i]->nucpos);
        }

        // currchr records current chromosome
        status = snprintf(currchr, sizeof currchr, "%s", r[0]->chr);
        if(status >= sizeof currchr) {
            fprintf(stderr, "%s:%d: buffer overflow\n", __FILE__, __LINE__);
            return BUFFER_OVERFLOW;
        }
        // Now get them all on the same position. Have to keep
        // checking chr in case one file moves to another chromosome.
        for(i = 0; onSameChr && i < n; ++i) {
            // Increment each reader so long as we're all on the same
            // chromosome and the reader's nucpos is low.
            while(onSameChr && r[i]->nucpos < maxnuc) {
                if( (status = DAFReader_next(r[i])) )
                    return status;
                diff = strcmp(r[i]->chr, currchr);
                if(diff != 0) {
                    // Assertion should succeed because DAFReader_next
                    // guarantees that chromosomes are in sort order.
                    assert(diff > 0);
                    onSameChr = 0;
                    imaxchr = i;
                }
            }
        }
    }
    while(!onSameChr || minnuc != maxnuc);

    // Make sure reference allele isn't fixed in readers. If it's
    // fixed, then we can't call the ancestral allele.
    double maxp, minp;
    maxp = minp = DAFReader_daf(r[0]);
    for(i=1; i < n; ++i) {
        double p = DAFReader_daf(r[i]); // derived allele freq
        minp = fmin(minp, p);
        maxp = fmax(maxp, p);
    }
    if(maxp == 0.0 || minp == 1.0)
        return NO_ANCESTRAL_ALLELE;

    if(!DAFReader_allelesMatch(n, r))
        return ALLELE_MISMATCH;

    return 0;
}

#endif

/// Return 1 if ancestral and derived alleles of all readers match; 0
/// otherwise
int DAFReader_allelesMatch(int n, DAFReader * r[n]) {
    char  *aa = r[0]->aa;
    char  *da = r[0]->da;
    int    daMissing = (0==strcmp(".", da));
    int    aaMissing = (0==strcmp(".", aa));
    int   i;
    for(i = 1; i < n; ++i) {
        if(0!=strcmp(aa, r[i]->aa))
            return 0;
        int curr_da_missing = (0 == strcmp(".", r[i]->da));
        int curr_aa_missing = (0 == strcmp(".", r[i]->aa));
        if(daMissing && !curr_da_missing) {
            daMissing = 0;
            da = r[i]->da;
            strcpy(r[0]->da, da);
        }
        if(aaMissing && !curr_aa_missing) {
            aaMissing = 0;
            aa = r[i]->aa;
            strcpy(r[0]->aa, aa);
        }
        if(!daMissing && !curr_da_missing
           && 0!=strcmp(da, r[i]->da))
            return 0;
        if(!aaMissing && !curr_aa_missing
           && 0!=strcmp(aa, r[i]->aa))
            return 0;
    }
    return 1;
}

/// Print header for daf file.
void DAFReader_printHdr(FILE * fp) {
    fprintf(fp, "%50s %5s %10s %2s %2s %8s\n",
            "file", "chr", "pos", "aa", "da", "daf");
}

/// Print current line of daf file
void DAFReader_print(DAFReader * r, FILE * fp) {
    assert(r->fname);
    fprintf(fp, "%50s %5s %10lu %2s %2s %8.6lf\n",
            r->fname, r->chr, r->nucpos, r->aa, r->da, r->p);
}

/// Return derived allele frequency of current line of daf file.
double DAFReader_daf(DAFReader * r) {
    assert(r->p >= 0.0 && r->p <= 1.0);
    return r->p;
}
