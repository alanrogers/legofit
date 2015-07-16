#include "gptree.h"

// One line of input for each segment of population tree.
// Each line is of form
//
//   x < 0.8 y < 0.2 z 2N=1e5 t=1000 samples=2
//
// which says that population x got 80% of its genes from y and 20%
// from z, has size 2N=1e5, duration 1000 generations, and that we
// sampled 2 haploid genomes.

int getDbl(double *x, Tokenizer *tkz, int i);

int getDbl(double *x, Tokenizer *tkz, int i) {
    char *end=NULL;

    *x = strtod(Tokenizer_token(tkz, i), &end);
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
    double p1, p2;
    int got_p1, got_p2, curr, nsamples;

    while(1) {
        curr = 0;
        got_p1 = got_p2 = nsamples = 0;
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
        name0 = Tokenizer_token(tkz, curr++);
        if('<' != *Tokenizer_token(tkz, curr++)) {
            fprintf(stderr,"Bad input line:\n");
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        if(!getDbl(&p1, tkz, curr)) {
            got_p1 = 1;
            if(p1 <= 0 || p1 >= 1.0) {
                fprintf(stderr,"Bad probability, %lf, in input line:\n", p1);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            ++curr;
        }
        name1 = Tokenizer_token(tkz, curr++);
        if('<' != *Tokenizer_token(tkz, curr++)) {
            fprintf(stderr,"Bad input line:\n");
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        if(!getDbl(&p2, tkz, curr)) {
            got_p2 = 1;
            if(p2 <= 0 || p2 >= 1.0) {
                fprintf(stderr,"Bad probability, p2=%lf, in input line:\n", p1);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            if(fabs(1.0 - p1 - p2) > 0.02) {
                fprintf(stderr,"%lf + %ld != 1.0\n", p1, p2);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            ++curr;
        }
        name2 = Tokenizer_token(tkz, curr++);
        if(0 != strcmp("2N", Tokenizer_token(tkz, curr++))) {
            fprintf(stderr,"Got %s when expecting 2N on input:\n",
                    Tokenizer_token(tkz, curr-1));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        if(!getDbl(&twoN, tkz, curr)) {
            if(twoN <= 0.0) {
                fprintf(stderr,"Bad value 2N=%lf, in input line:\n", twoN);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            ++curr;
        }
    }

    Tokenizer_free(tkz);

    return root;
}
