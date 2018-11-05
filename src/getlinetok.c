#include "getlinetok.h"
#include "misc.h"
#include "tokenizer.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct GetLineTok {
    size_t buffsize;
    char *buff;
    Tokenizer *tkz;
    FILE *fp; // not locally owned
};

static int strempty(const char *s);

/// Return true if string contains only whitespace; false otherwise.
static int strempty(const char *s) {
    const char *p = s;

    while(isspace(*p))
        ++p;
    if(*p == '\0')
        return true;
    return false;
}

GetLineTok *GetLineTok_new(size_t buffsize, int maxtokens, FILE *fp) {
    if(buffsize < 2)
        DIE("buffsize argument must be >1");
    GetLineTok *self = malloc(sizeof(GetLineTok));
    if(self==NULL)
        DIE("bad malloc");
    self->buffsize = buffsize;
    self->buff = malloc(buffsize);
    if(self->buff == NULL)
        DIE("bad malloc");
    self->buff[0] = '\0';
    self->fp = fp;
    self->tkz = Tokenizer_new(maxtokens);
    return self;
};

/// Return 0 on success, EOF or ENOMEM on failure.
int GetLineTok_next(GetLineTok *self, const char *sep,
                        const char *extraneous) {
    self->buff[0] = '\0';
    Tokenizer_clear(self->tkz);
    do {
        if(NULL == fgets(self->buff, self->buffsize, self->fp))
            return EOF;
    } while(strempty(self->buff));

    while(!strchr(self->buff, '\n') && !feof(self->fp)) {
        assert(strlen(self->buff) == self->buffsize-1);
        // buffer overflow: reallocate
        size_t oldsize = strlen(self->buff);
        self->buffsize *= 2;
        char *oldbuff = self->buff;
        self->buff = realloc(self->buff, self->buffsize);
        if(self->buff == NULL) {
            free(oldbuff);
            return ENOMEM;
        }
        // Try to read the rest of the current line. If there is
        // nothing more to read, that's OK, because we've already
        // filled the first half of self->buff.
        if(NULL == fgets(self->buff + oldsize,
                         self->buffsize - oldsize,
                         self->fp))
            break;
    }

    Tokenizer_split(self->tkz, self->buff, sep);
    Tokenizer_strip(self->tkz, extraneous);

    return 0;
}

int GetLineTok_ntokens(GetLineTok *self) {
    return Tokenizer_ntokens(self->tkz);
}

char *GetLineTok_token(GetLineTok * self, int index) {
    return Tokenizer_token(self->tkz, index);
}

void GetLineTok_print(GetLineTok *self, FILE *ofp) {
    Tokenizer_print(self->tkz, ofp);
}
