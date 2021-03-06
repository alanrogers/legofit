/**
 * @file tokenizer.c
 * @author Alan R. Rogers
 * @brief Tokenize a character string.
 *
 * This file implements a class that tokenizes a string. Usage is like
 * this:
 *
 *    int ntokens;
 *    char buff[100];
 *    Tokenizer *tkz = Tokenizer_new(maxTokens);
 *
 *    strcpy(buff, "my: input   ; ; string");
 *    ntokens = Tokenizer_split(tkz, buff, ":;");
 *
 * The 3rd argument to Tokenizer_split defines the characters
 * that separate tokens. Then you can access individual tokens like
 * this:
 *
 *    char *token;
 *
 *    token = Tokenizer_token(tkz, 3);
 *
 * The tokens themselves still reside in buff, so this array must not
 * change until the next call to Tokenizer_split.
 *
 * The memory allocated by Tokenizer_new is freed by Tokenizer_free.
 *
 * The argument to Tokenizer_new determines the initial size of an
 * internal array of pointers to tokens. If Tokenizer_split is handed
 * a string with more than this number of tokens, it will re-allocate
 * the internal array. It is legal to initialize with
 * Tokenizer_new(0).
 *
 * @copyright Copyright (c) 2014,2018 Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "misc.h"
#include "tokenizer.h"

/// Tokenizer constructor
Tokenizer  *Tokenizer_new(int maxTokens) {
    Tokenizer  *t = malloc(sizeof(Tokenizer));

    CHECKMEM(t);

    if(maxTokens == 0)
        t->tokptr = NULL;
    else {
        t->tokptr = malloc(maxTokens * sizeof(t->tokptr[0]));
        CHECKMEM(t->tokptr);
    }

    t->maxTokens = maxTokens;
    return t;
}

/// Tokenizer destructor
void Tokenizer_free(Tokenizer * self) {
    assert(self);
    if(self->tokptr)
        free(self->tokptr);
    free(self);
}

/**
 * Turn string "buff" into an array of tokens, assuming that tokens in
 * the input string may be separated by any of the characters in
 * string "sep". Supress empty tokens. Return the number of tokens.
 */
int Tokenizer_split(Tokenizer * self, char *buff, const char *sep) {
    assert(self);
    assert(buff);
    assert(sep);
    char *ptr = buff, *token;

    self->n = 0;
    while((token = strsep(&ptr, sep)) != NULL) {
        if(*token == '\0')
            continue;           // skip empty tokens
        if(self->n == self->maxTokens) {
            // reallocate
            int need = 1;
            if(ptr)
                need += strCountSetChunks(ptr, sep);
            self->maxTokens += need;
            self->tokptr = realloc(self->tokptr,
                                   self->maxTokens * sizeof(self->tokptr[0]));
            if(self->tokptr == NULL)
                DIE("bad realloc");
        }
        assert(self->n < self->maxTokens);
        self->tokptr[self->n] = token;
        ++(self->n);
    }

    return self->n;
}

/// Return pointer to token with given index
char       *Tokenizer_token(Tokenizer * t, int ndx) {
    assert(t);
    if(t->n == 0)
        DIE("tokenizer has 0 tokens");
    if(ndx >= t->n) {
        fprintf(stderr, "%s:%s:%d: ndx=%d, but there are only %d tokens.\n",
                __FILE__, __func__, __LINE__, ndx, t->n);
        exit(EXIT_FAILURE);
    }
    return t->tokptr[ndx];
}

/**
 * Strip extraneous chars (those list in "extraneous") from both ends
 * of each token. If the result is an empty string, this token is
 * removed from the list. The function returns the number of tokens.
 */
int Tokenizer_strip(Tokenizer * t, const char *extraneous) {
    assert(t);
    int         curr = 0;
    char       *end;

    while(curr < t->n) {

        /* move beginning beyond initial extraneous chars */
        t->tokptr[curr] += strspn(t->tokptr[curr], extraneous);

        /* set end = terminating '\0' of current token */
        end = rindex(t->tokptr[curr], '\0');

        while(t->tokptr[curr] < end) {
            /* Is last char extraneous? */
            if(strchr(extraneous, *(end - 1)) == NULL) {
                /* Last char not extraneous. We're done. */
                break;
            } else {
                /* last char in token is extraneous: remove it */
                --end;
                *end = '\0';
            }
        }

        if(t->tokptr[curr] == end) {
            /* remove empty token from list */
            assert(strlen(t->tokptr[curr]) == 0);
            memmove(t->tokptr + curr, t->tokptr + curr + 1,
                    (t->n - curr - 1) * sizeof(t->tokptr[0]));
            --(t->n);
        } else
            ++curr;
    }
    return t->n;
}

/// Return number of tokens
int Tokenizer_ntokens(Tokenizer * t) {
    assert(t);
    return t->n;
}

/**
 * Search for string s among tokens. On success, return index of token.
 * On failure, return the current number of tokens. After each call,
 * the returned value should be compared with that of
 * Tokenizer_ntokens.
 */
int Tokenizer_find(Tokenizer * t, const char *s) {
    assert(t);
    assert(s);
    int         i;

    for(i = 0; i < t->n; ++i) {
        if(strcmp(t->tokptr[i], s) == 0)
            return i;
    }
    return t->n;
}

/// Print a summary of the information in a Tokenizer
void Tokenizer_printSummary(const Tokenizer * tkz, FILE * ofp) {
    assert(tkz);
    assert(ofp);
    fprintf(ofp,
            "Tokenizer:         n=%8d (number of tokens)\n"
            "           maxTokens=%8d\n", tkz->n, tkz->maxTokens);
}

/// Print Tokenizer object
void Tokenizer_print(const Tokenizer * tkz, FILE * ofp) {
    assert(tkz);
    assert(ofp);
    int         i;

    for(i = 0; i < tkz->n; ++i)
        fprintf(ofp, " \"%s\"", tkz->tokptr[i]);
    putc('\n', ofp);
}

void Tokenizer_clear(Tokenizer *self) {
    self->n = 0;
}
