#include "gptree.h"
#include "tokenizer.h"

// One line of input for each segment of population tree.
// Each line is of form
//
//   x < 0.8 y < z 2N=1e5 t=1000 samples=2
//
// which says that population x got 80% of its genes from y and 20%
// from z, has size 2N=1e5, duration 1000 generations, and that we
// sampled 2 haploid genomes. If there was only a single parental
// population, the input line would look like this
//
//   x < y 2N=1e5 t=1000 samples=2

#define CHECK_INDEX(ndx,n) do{                                          \
    if((ndx)>=(n)){                                                     \
    fprintf(stderr,"%s:%s:%d: index out of bounds\n",                   \
            __FILE__,__func__,__LINE__);                                \
    exit(EXIT_FAILURE);                                                 \
    }while(0)

int getDbl(double *x, Tokenizer *tkz, int i);
int getULong(unsigned long *x, Tokenizer *tkz, int i);

int getDbl(double *x, Tokenizer *tkz, int i) {
    char *end=NULL;

    *x = strtod(Tokenizer_token(tkz, i), &end);
    if(end!=NULL && *end == '\0')
        return 0;  // success
    return 1;      // failure
}

int getULong(unsigned long *x, Tokenizer *tkz, int i) {
    char *end=NULL;

    *x = strtoul(Tokenizer_token(tkz, i), &end, 10);
    if(end!=NULL && *end == '\0')
        return 0;  // success
    return 1;      // failure
}


PopNode *mktree(FILE *fp, HashTab *ht) {
    int ntokens;
    char buff[500];
    Tokenizer *tkz = Tokenizer_new(50);

    char *name0, *name1, *name2, *s;
    El *el0, *el1, *el2;
    double p, twoN, t;
    unsigned samples;
    int got_p, curr, nsamples;

    while(1) {
        curr = 0;
        got_p = nsamples = 0;
        name0 = name1 = name2 = NULL;
        el0 = el1 = el2 = NULL;
        if(fgets(buff, sizeof buff, fp) == NULL)
            break;

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n",
                    __FILE__, __LINE__, sizeof(buff));
        Tokenizer_split(tkz, buff, " \t="); // tokenize
        ntokens = Tokenizer_strip(tkz, " \t=\n");
        CHECK_INDEX(curr, ntokens);
        name0 = Tokenizer_token(tkz, curr++);
        CHECK_INDEX(curr, ntokens);
        if('<' != *Tokenizer_token(tkz, curr++)) {
            fprintf(stderr,"Bad input line:\n");
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        CHECK_INDEX(curr, ntokens);
        if(!getDbl(&p, tkz, curr)) {
            got_p = 1;
            if(p <= 0 || p >= 1.0) {
                fprintf(stderr,"Bad probability, %lf, in input line:\n", p1);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            ++curr;
        }
        CHECK_INDEX(curr, ntokens);
        name1 = Tokenizer_token(tkz, curr++);
        CHECK_INDEX(curr, ntokens);
        if('<' == *Tokenizer_token(tkz, curr)) {
            name2 = Tokenizer_token(tkz, ++curr);
            ++curr;
        }
        CHECK_INDEX(curr, ntokens);
        if(0 != strcmp("2N", Tokenizer_token(tkz, curr++))) {
            fprintf(stderr,"Got %s when expecting 2N on input:\n",
                    Tokenizer_token(tkz, curr-1));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        CHECK_INDEX(curr, ntokens);
        if(!getDbl(&twoN, tkz, curr)) {
            if(twoN <= 0.0) {
                fprintf(stderr,"Bad value 2N=%lf, in input line:\n", twoN);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            ++curr;
        }else{
            fprintf(stderr,"Can't parse \"%s\" as a double. Expecting value of 2N\n",
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        CHECK_INDEX(curr, ntokens);
        if(0 != strcmp("t", Tokenizer_token(tkz, curr++))) {
            fprintf(stderr,"Got %s when expecting t on input:\n",
                    Tokenizer_token(tkz, curr-1));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        CHECK_INDEX(curr, ntokens);
        if(!getDbl(&t, tkz, curr)) {
            if(t <= 0.0) {
                fprintf(stderr,"Bad value t=%lf, in input line:\n", t);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            ++curr;
        }else{
            fprintf(stderr,"Can't parse \"%s\" as a double. Expecting value of t\n",
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        if(curr == ntokens)
            continue;
        CHECK_INDEX(curr, ntokens);
        if(0 != strcmp("samples", Tokenizer_token(tkz, curr++))) {
            fprintf(stderr,"Got %s when expecting \"samples\" on input:\n",
                    Tokenizer_token(tkz, curr-1));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        CHECK_INDEX(curr, ntokens);
        if(!getULong(&samples, tkz, curr)) {
            if(samples > 1000) {
                fprintf(stderr,"Bad value samples=%u, in input line:\n", samples);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            ++curr;
        }else{
            fprintf(stderr,"Can't parse \"%s\" as an unsigned int."
                    " Expecting value of \"samples\"\n",
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        if(curr < ntokens) {
            fprintf(stderr,"Too many tokens on input line:\n");
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
    }

    Tokenizer_free(tkz);

    return root;
}
