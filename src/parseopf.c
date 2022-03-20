/**
 * @file parseopf.c
 * @author Alan R. Rogers
 * @brief Parse a .opf file, containing site pattern frequencies.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "parseopf.h"
#include "branchtab.h"
#include "lblndx.h"
#include "misc.h"
#include "string.h"
#include "tokenizer.h"
#include <stdio.h>

/// Parse a .opf file containing observed site pattern frequencies.
/// Recognizes comments, which extend from '#' to end-of-line.
BranchTab *parseOpf(const char *fname, const LblNdx *lblndx) {
    FILE *fp = efopen(fname, "r");

    BranchTab *self = BranchTab_new(LblNdx_size(lblndx));
    CHECKMEM(self);

    int         ntokens;
    char        buff[500];
    Tokenizer  *tkz = Tokenizer_new(50);

    while(1) {
        if(fgets(buff, sizeof(buff), fp) == NULL)
            break;

        if(!strchr(buff, '\n') && !feof(fp))
            eprintf("s:%s:%d: buffer overflow. buff size: %zu\n",
                    __FILE__, __func__, __LINE__, sizeof(buff));

        // strip trailing comments
        char *comment = strchr(buff, '#');
        if(comment)
            *comment = '\0';

        Tokenizer_split(tkz, buff, " \t"); // tokenize
        ntokens = Tokenizer_strip(tkz, " \t\n");
        if(ntokens == 0)
            continue;

        char *tok = Tokenizer_token(tkz, 0);
        tipId_t key=LblNdx_getTipId(lblndx, tok);
        if(key==0) {
            fprintf(stderr,"%s:%s:%d: can't find id for label %s\n",
                    __FILE__,__func__,__LINE__, tok);
            fprintf(stderr,"  Either there is no \"segment v\" in the"
                    " .lgo file, or that segment has no samples.\n");
            exit(EXIT_FAILURE);
        }

        tok = Tokenizer_token(tkz, 1);
        errno = 0;
        double prob = strtod(tok, NULL);
        if(errno) {
            fprintf(stderr,"%s:%s:%d: Can't parse 2nd field as float.\n",
                    __FILE__,__func__,__LINE__);
            fprintf(stderr," input:");
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }

        BranchTab_add(self, key, prob);
    }
    fclose(fp);
    Tokenizer_free(tkz);
    return self;
}
