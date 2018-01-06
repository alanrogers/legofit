/**
@file scrmreader.c
@page scrmreader
@brief Interface to scrm output files.
*/

#include "scrmreader.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <errno.h>

struct ScrmReader {
    int npops;
    int *nsamples;
    double *daf;

    // Independent replicates in scrm output appear as separate
    // chromosomes, which are labelled by unsigned integers.
    unsigned chr;

    // Position values in scrm output are ignored. Instead, ScrmReader
    // returns positions as a sequence of unsigned integers.
    unsigned long nucpos;
    FILE *fp;
};

// Allocate and initialize a new ScrmReader from input stream.
ScrmReader *ScrmReader_new(FILE *fp) {
    // buffer is large, because scrm command lines can be long
    char buff[8192];
    errno = 0;
    if(fgets(buff, sizeof buff, fp) == NULL) 
        return NULL;

    if(!strchr(buff, '\n') && !feof(fp)) {
        fprintf(stderr, "%s:%d: buffer overflow. buff size: %zu\n",
                __FILE__, __LINE__, sizeof buff);
        fprintf(stderr,"input: %s\n", buff);
        return NULL;
    }

    ScrmReader *self = malloc(sizeof ScrmReader);
    CHECKMEM(self);
    memset(self, 0, sizeof ScrmReader);
    self->fp = fp;

    Tokenizer *tkz = Tokenizer_new((sizeof buff)/2);
    Tokenizer_split(tkz, buff, " ");
    int ntokens = Tokenizer_strip(tkz, " \n");

    if(strcmp("scrm", Tokenizer_token(tkz, 0)) != 0) {
        fprintf(stderr,"%s:%d: input file is not scrm output\n",
                __FILE__,__LINE__);
        free(self);
        exit(EXIT_FAILURE);
    }

    int i, j, npops;
    long h;
    char *token, *end;

    // 
    for(i=1; i < ntokens; ++i) {
        token = Tokenizer_token(tkz, i);
        if(strcmp("-I", token) == 0) {
            token = Tokenizer_token(tkz, i+1);
            if(self->npops == 0) {
                // set npops and allocate nsamples
                self->npops = strtol(token, NULL, 10);
                self->nsamples = malloc(self->npops * sizeof(int));
                CHECKMEM(self->nsamples);
                for(j=0; j < self->npops; ++j) {
                    token = Tokenizer_token(tkz, i+2+j);
                    self->nsamples[j] = strtol(token, NULL, 10);
                }
            }else{
                // check for consistency and increment nsamples
                npops = strtol(token, NULL, 10);
                if(npops != self->npops) {
                    fprintf(stderr,"%s:%d: inconsistent population count:"
                            " %d != %d\n",
                            __FILE__,__LINE__,npops, self->npops);
                    exit(EXIT_FAILURE);
                }
                for(j=0; j < self->npops; ++j) {
                    token = Tokenizer_token(tkz, i+2+j);
                    self->nsamples[j] += strtol(token, NULL, 10);
                }
            }
            i += self->npops;
        }else if(strcmp("-eI", token) == 0) {
            if(self->npops == 0) {
                // count populations and allocate nsamples
                for(j=i+2; j<ntokens; ++j) {
                    token = Tokenizer_token(tkz, j);
                    h = strtol(token, &end, 10);
                    if(end==token) // token isn't an integer
                        break;
                    else           // token is an integer
                        ++self->npops;
                }
                assert(self->nsamples == NULL);
                if(self->npops == 0) {
                    fprintf(stderr,"%s:%d: npops is zero\n",
                            __FILE__,__LINE__);
                    exit(EXIT_FAILURE);
                }
                self->nsamples = malloc(self->npops * sizeof(int));
                for(j=0; j<self->npops; ++j) {
                    token = Tokenizer_token(tkz, i+2+j);
                    h = strtol(token, &end, 10);
                    assert(end != token);
                    self->nsamples[j] = h;
                }
            }else{
                // increment nsamples
                assert(self->npops > 0);
                assert(self->nsamples == NULL);
                for(j=0; j<self->npops; ++j) {
                    token = Tokenizer_token(tkz, i+2+j);
                    h = strtol(token, &end, 10);
                    assert(end != token);
                    self->nsamples[j] += h;
                }
            }
        }
    }
    
}

// Rewind input and reset chr and nucpos. Doesn't work
// if input is stdin.
int ScrmReader_rewind(ScrmReader *self);

// Move ScrmReader to next nucleotide site.
int ScrmReader_next(ScrmReader *self);

// Return current chromosome.
unsigned ScrmReader_chr(ScrmReader *self);

// Return current nucleotide position.
unsigned long ScrmReader_nucpos(ScrmReader *self);

// Return number of populations.
int ScrmReader_npops(ScrmReader *self);

// Return number of samples from population i.
int ScrmReader_nsamples(ScrmReader *self, int i);

// Return frequency of derived allele in sample from population i.
double ScrmReader_daf(ScrmReader *self, int i);
