/**
@file simreader.c
@page simreader
@brief Interface to output files produced by ms2sim and msprime.

There are two versions of the input format. In version 1, the top of
the input looked like this:

    npops = 4
    pop sampsize
    x 1
    y 1
    n 1
    d 1
    0 0 0 1 1 
    0 0 0 1 0

In version 2, it looks like this:

    nchromosomes = 10
    chrlen = 2000000
    npops = 5
    pop sampsize
    x 1
    y 1
    n 1
    d 1
    s 1
    0 345.7722022796852 0 0 0 1 0 
    0 377.99065190067955 1 0 0 0 0 
    0 407.9258238758751 0 0 0 0 1 

Version 2 has 2 lines at the top, giving the number of chromosomes and
the length in basepairs of each chromosome. In the lines of data that
follow the header, the 2nd field has the nucleotide position of the
current site.

The SimReader_new function distinguishes between these formats by
examining the first character string in the file. If that string is
"npops", we are reading format 1. If it's "nchromosomes", we're in
format 2.

If the input data are in format 2, the function SimReader_pos returns
the position of the current site. If the input is format 1, this
function returns -1.
*/

#include "simreader.h"
#include "misc.h"
#include "tokenizer.h"
#include "error.h"
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>

struct SimReader {
    int input_format; // 1 or 2
    int nchromosomes;
    int chr_len; // length of chromosome in basepairs.
    int sampleDim;
    unsigned *nsamples;
    char **lbl;
    double *daf;
    Tokenizer *tkz;
    FILE *fp;  // not locally owned

    // Independent replicates in output appear as separate
    // chromosomes, which are labelled by integers.
    int chr;

    // Nucleotide position within current chromosome
    double pos;

    // number of fields in input data
    int nfields;
};

// destructor doesn't free self->fp
void SimReader_free(SimReader *self) {
    assert(self);
    free(self->nsamples);
    for(int i=0; i < self->sampleDim; ++i)
        free(self->lbl[i]);
    free(self->lbl);
    free(self->daf);
    Tokenizer_free(self->tkz);
    free(self);
}

// Allocate and initialize a new SimReader from input stream.
SimReader *SimReader_new(FILE *fp) {

    char buff[128];
    int i, status;

    SimReader *self = malloc(sizeof(SimReader));
    CHECKMEM(self);
    memset(self, 0, sizeof(SimReader));

    self->chr = -1;
    self->pos = -1.0;
    self->fp = fp;
    assert(fp != NULL);

    self->tkz = Tokenizer_new(sizeof(buff)/2);

    status = readline(sizeof(buff), buff, self->fp);
    switch(status) {
    case 0:
        break;
    case EOF:
        fprintf(stderr,"%s:%d: unexpected EOF\n", __FILE__,__LINE__);
        goto new_fail;
    case BUFFER_OVERFLOW:
        fprintf(stderr,"%s:%d: buffer overflow reading msprime command\n",
                __FILE__,__LINE__);
        goto new_fail;
    default:
        fprintf(stderr,"%s:%d: unknown error\n", __FILE__,__LINE__);
        goto new_fail;
    }

    int ntokens;

    // parse header
    Tokenizer_split(self->tkz, buff, "=");
    ntokens = Tokenizer_strip(self->tkz, " \n");
    if(ntokens != 2)
        goto bad_header;

    if( 0 == strcmp("nchromosomes", Tokenizer_token(self->tkz, 0)) ) {
        self->input_format = 2;

        // set nchromosomes
        self->nchromosomes = strtol(Tokenizer_token(self->tkz, 1), NULL, 10);
        if(self->nchromosomes <= 0)
            goto bad_header; 

        // set chr_len
        status = readline(sizeof(buff), buff, self->fp);
        if(status)
            goto bad_header;
        Tokenizer_split(self->tkz, buff, "=");
        ntokens = Tokenizer_strip(self->tkz, " \n");
        if(ntokens != 2)
            goto bad_header;
        if( 0 != strcmp("chrlen", Tokenizer_token(self->tkz, 0)) )
            goto bad_header;
        self->chr_len = strtol(Tokenizer_token(self->tkz, 1), NULL, 10);
        if(self->chr_len <= 0)
            goto bad_header; 

        // set sampleDim
        status = readline(sizeof(buff), buff, self->fp);
        if(status)
            goto bad_header;
        Tokenizer_split(self->tkz, buff, "=");
        ntokens = Tokenizer_strip(self->tkz, " \n");
        if(ntokens != 2)
            goto bad_header;
        if( 0 != strcmp("npops", Tokenizer_token(self->tkz, 0)) )
            goto bad_header;
        self->sampleDim = strtol(Tokenizer_token(self->tkz, 1), NULL, 10);
        if(self->sampleDim <= 0)
            goto bad_header; 
        
    }else if( 0 == strcmp("npops", Tokenizer_token(self->tkz, 0)) ) {
        self->input_format = 1;

        // set sampleDim
        self->sampleDim = strtol(Tokenizer_token(self->tkz, 1), NULL, 10);
        if(self->sampleDim <= 0)
            goto bad_header;
    }else
        goto bad_header;

    self->nsamples = malloc(self->sampleDim * sizeof(self->nsamples[0]));
    CHECKMEM(self->nsamples);
    self->lbl = malloc(self->sampleDim * sizeof(self->lbl[0]));
    CHECKMEM(self->lbl);
    memset(self->lbl, 0, self->sampleDim * sizeof(self->lbl[0]));

    // read population labels and sample sizes
    status = readline(sizeof(buff), buff, self->fp);
    if(status)
        goto new_fail;
    Tokenizer_split(self->tkz, buff, " ");
    ntokens = Tokenizer_strip(self->tkz, " \n");
    if(ntokens != 2
       || 0!=strcmp("pop", Tokenizer_token(self->tkz, 0))
       || 0!=strcmp("sampsize", Tokenizer_token(self->tkz, 1))){
        fprintf(stderr,"%s:%d: error reading header\n",
                __FILE__,__LINE__);
        goto new_fail;
    }

    // Label and size of sample for each population
    for(i=0; i < self->sampleDim; ++i) {
        status = readline(sizeof(buff), buff, self->fp);
        switch(status) {
        case 0:
            break;
        case EOF:
            fprintf(stderr,"%s:%d: unexpected EOF\n", __FILE__,__LINE__);
            goto new_fail;
        case BUFFER_OVERFLOW:
            fprintf(stderr,"%s:%d: buffer overflow\n", __FILE__,__LINE__);
            goto new_fail;
        default:
            fprintf(stderr,"%s:%d: unknown error\n", __FILE__,__LINE__);
            goto new_fail;
        }
        Tokenizer_split(self->tkz, buff, " ");
        ntokens = Tokenizer_strip(self->tkz, " \n");
        if(ntokens != 2) {
            fprintf(stderr,"%s:%d: expecting 2 tokens; got %d\n",
                    __FILE__,__LINE__, ntokens);
            exit(EXIT_FAILURE);
        }

        self->lbl[i] = strdup(Tokenizer_token(self->tkz, 0));
        self->nsamples[i] = strtol(Tokenizer_token(self->tkz, 1),
                                   NULL, 10);
    }

    // allocate daf array
    self->daf = malloc(self->sampleDim * sizeof(self->daf[0]));
    CHECKMEM(self->daf);

    // number of fields in each line of input data.
    self->nfields = 1; // first field is chromosome
    if(self->input_format == 2)
        self->nfields += 1; // 2nd field is position in format 2
    for(i=0; i < self->sampleDim; ++i)  // others are haploid genotypes
        self->nfields += self->nsamples[i];

    return self;

 bad_header:    
    fprintf(stderr,"%s:%d: bad input header\n",
            __FILE__,__LINE__);
    Tokenizer_print(self->tkz, stderr);
    if(self->tkz)
        Tokenizer_free(self->tkz);
    free(self);
    return NULL;

 new_fail:
    if(self==NULL)
        return NULL;
    if(self->nsamples)
        free(self->nsamples);
    if(self->lbl) {
        for(i=0; i < self->sampleDim; ++i) {
            if(self->lbl[i])
                free(self->lbl[i]);
        }
        free(self->nsamples);
    }
    if(self->daf)
        free(self->daf);
    if(self->tkz)
        Tokenizer_free(self->tkz);
    free(self);
    return NULL;
}

// Move SimReader to next nucleotide site.
int SimReader_next(SimReader *self) {
    char buff[128];
    int status, ntokens=0;
    do{
        status = readline(sizeof(buff), buff, self->fp);
        if(status)
            return status;
        Tokenizer_split(self->tkz, buff, " ");
        ntokens = Tokenizer_strip(self->tkz, " \n");
    }while(ntokens==0 && !feof(self->fp));

    if(ntokens != self->nfields) {
        fprintf(stderr,"%s:%d: read %d tokens; was expecting %d\n",
                __FILE__,__LINE__, ntokens, self->nfields);
        exit(EXIT_FAILURE);
    }

    // current chromosome
    self->chr = strtol(Tokenizer_token(self->tkz, 0), NULL, 10);

    int start; // index of first field of genotypes

    // set start index of genotypes and (format 2) nucleotide position
    if(self->input_format == 1) {
        start=1;   // skip chr, in position 0
    }else{
        assert(self->input_format == 2);
        start=2;   // skip chr and pos
        self->pos = strtod(Tokenizer_token(self->tkz, 1), NULL);
    }

    // calculate derived allele frequency w/i each pop
    double nderived;
    int pop, i;
    char *token, *end;
    for(pop=0; pop < self->sampleDim; ++pop) {
        nderived = 0.0;
        for(i=start; i < start + self->nsamples[pop]; ++i) {
            if(i >= Tokenizer_ntokens(self->tkz)) {
                fprintf(stderr,"%s:%s:%d: too few genotypes\n",
                        __FILE__,__func__,__LINE__);
                return EDOM;
            }
            token = Tokenizer_token(self->tkz, i);
            unsigned gtype = strtoul(token, &end, 10);
            if(token==end || (gtype!=0 && gtype!=1)) {
                fprintf(stderr,"%s:%d: illegal genotype: %s\n",
                        __FILE__,__LINE__, token);
                return EDOM;
            }
            nderived += gtype;
        }
        self->daf[pop] = nderived / self->nsamples[pop];
        start += self->nsamples[pop];
    }
    return 0;
}

// Return current chromosome.
unsigned SimReader_chr(const SimReader *self) {
    return self->chr;
}

// Return current nucleotide position if input data are in format 2.
// Return -1.0 if input is format 1.
double SimReader_pos(const SimReader *self) {
    return self->pos;
}

// Return number of populations.
int SimReader_sampleDim(SimReader *self) {
    return self->sampleDim;
}

// Return number of samples from population i.
int SimReader_nsamples(SimReader *self, int i) {
    assert(i < self->sampleDim);
    return self->nsamples[i];
}

// Return frequency of derived allele in sample from population i.
double SimReader_daf(SimReader *self, int i) {
    assert(i < self->sampleDim);
    assert(self->chr >= 0);
    return self->daf[i];
}

const char *SimReader_lbl(SimReader *self, int i) {
    assert(i < self->sampleDim);
    assert(self->chr >= 0);
    return self->lbl[i];
}

#ifdef TEST

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

const char *tstInput =
    "npops = 4\n"
    "pop sampsize\n"
    "x 1\n"
    "y 1\n"
    "n 1\n"
    "d 1\n"
    "0 1 0 0 0 \n"
    "0 0 1 0 0 \n"
    "1 0 1 1 1 \n"
    "1 1 0 0 0 \n";

int main(int argc, char **argv) {
    int         verbose = 0;
    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xmspreader [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    const char *tstfname = "mspreader.tst";
    FILE *fp = fopen(tstfname, "w");
    if(fp==NULL)
        DIE("bad fopen");
    fputs(tstInput, fp);
    fclose(fp);

    fp = fopen(tstfname, "r");
    SimReader *r = SimReader_new(fp);
    CHECKMEM(r);

    assert(-1 == SimReader_chr(r));
    int dim = SimReader_sampleDim(r);
    assert(4 == dim);

    int i, status;

    for(i=0; i<dim; ++i) {
        assert(1 == SimReader_nsamples(r, i));
    }

    // test mspreader

    status = SimReader_next(r);
    assert(status==0);
    assert(0 == SimReader_chr(r));
    assert(1.0 == SimReader_daf(r, 0));
    assert(0.0 == SimReader_daf(r, 1));
    assert(0.0 == SimReader_daf(r, 2));
    assert(0.0 == SimReader_daf(r, 3));
    if(verbose) {
        for(i=0; i<SimReader_sampleDim(r); ++i) {
            printf("lbl[%d]=%s ", i, SimReader_lbl(r, i));
            printf("daf[%d]=%g\n", i, SimReader_daf(r, i));
        }
    }

    status = SimReader_next(r);
    if(status) {
        fprintf(stderr,"%s:%d: ERR SimReader_next returned %d\n",
                __FILE__,__LINE__, status);
        exit(EXIT_FAILURE);
    }
    assert(0 == SimReader_chr(r));
    assert(0.0 == SimReader_daf(r, 0));
    assert(1.0 == SimReader_daf(r, 1));
    assert(0.0 == SimReader_daf(r, 2));
    assert(0.0 == SimReader_daf(r, 3));

    status = SimReader_next(r);
    if(status) {
        fprintf(stderr,"%s:%d: ERR SimReader_next returned %d\n",
                __FILE__,__LINE__, status);
        exit(EXIT_FAILURE);
    }
    assert(1 == SimReader_chr(r));
    assert(0.0 == SimReader_daf(r, 0));
    assert(1.0 == SimReader_daf(r, 1));
    assert(1.0 == SimReader_daf(r, 2));
    assert(1.0 == SimReader_daf(r, 3));

    status = SimReader_next(r);
    if(status) {
        fprintf(stderr,"%s:%d: ERR SimReader_next returned %d\n",
                __FILE__,__LINE__, status);
        exit(EXIT_FAILURE);
    }
    assert(1 == SimReader_chr(r));
    assert(1.0 == SimReader_daf(r, 0));
    assert(0.0 == SimReader_daf(r, 1));
    assert(0.0 == SimReader_daf(r, 2));
    assert(0.0 == SimReader_daf(r, 3));

    fclose(fp);
    SimReader_free(r);
    remove(tstfname);
    unitTstResult("SimReader", "OK");

    return 0;
}

#endif
