/**
 * @file tokenizer.h
 * @author Alan R. Rogers
 * @brief Header for tokenizer.c
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#ifndef HAVE_TOKENIZER_H
#define HAVE_TOKENIZER_H
#include <stdio.h>
#include "typedefs.h"

/*
 * Defined here rather than hidden in tokenizer.c so that tokenizers
 * can be allocated as automatic variables, i.e. on the stack, rather
 * than via a call to malloc. This usage looks like:
 *
 *    Tokenizer tkz;
 *
 *    Tokenizer_split(&tkz, buff, sep);
 *
 * On a 32 bit machine, the max values of ints and longs are the same:
 * 2.14e9.  On a 64 bit machine, ints are the same, but longs are
 * 9.22e18. I'm using ints here, because I don't think I'll need to
 * tokenize more than 2.14e9 items.
 */
struct Tokenizer {
    int         n;              /* number of tokens */
    int         maxTokens;      /* maximum number of tokens */
    char      **tokptr;         /* array of pointers to tokens */
};

Tokenizer  *Tokenizer_new(int maxtok);
void        Tokenizer_free(Tokenizer * t);
int         Tokenizer_split(Tokenizer * t, char *buff, const char *sep);
char       *Tokenizer_token(Tokenizer * t, int index);
int         Tokenizer_ntokens(Tokenizer * t);
int         Tokenizer_strip(Tokenizer * t, const char *extraneous);
int         Tokenizer_find(Tokenizer * t, const char *s);
void        Tokenizer_print(const Tokenizer * tkz, FILE * ofp);
void        Tokenizer_printSummary(const Tokenizer * tkz, FILE * ofp);
#endif
