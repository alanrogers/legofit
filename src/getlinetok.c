#include "getlinetoc.h"
#include <stdio.h>
#include <stdlib.h>

struct GetLineTok {
    size_t buffsize;
    int maxtokens;
    char *buff;
    Tokenizer *tkz;
    FILE *fp; // not locally owned
};

int strcount(char *s, const char *sep);
int strempty(const char *s);

/// Return true if string contains only whitespace; false otherwise.
int strempty(const char *s) {
    const char *p = s;

    while(isspace(*p))
        ++p;
    if(*p == '\0')
        return true;
    return false;
}

/// Count the characters in string "s" that are members of string
/// "sep".
int strcount(char *s, const char *sep) {
    count=0;
    while(*s) {
        if(strchr(sep, *s))
            ++count;
        ++s;
    }
    return count;
}

GetLineTok *GetLineTok_new(size_t buffsize, int maxtokens, FILE *fp) {
    GetLineTok *self = malloc(sizeof(GetLineTok));
    if(self==NULL)
        DIE("bad malloc");
    new->buffsize = buffsize;
    self->buff = malloc(buffsize);
    if(self->buff == NULL)
        DIE("bad malloc");
    self->buff[0] = '\0';
    self->fp = fp;
    self->maxtokens = maxtokens;
    self->tkz = Tokenizer_new(maxtokens);
    return self;
};

/// Read a line and tokenize, reallocating if necessary.
/// Return 0 on success, EOF or ENOMEM on failure.
int GetLineTok_nextLine(GetLineTok *self, const char *sep,
                        const char *extraneous) {
    do {
        if(NULL == fgets(self->buff, self->buffsize, self->fp))
            return EOF;
    } while(strempty(buff));

    while(!strchr(self->buff, '\n') && !feof(self->ifp)) {
        assert(strlen(self->buff) == self->buffsize);
        // buffer overflow: reallocate
        size_t oldsize = self->buffsize;
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
        if(NULL == fgets(self->buff + oldsize, oldsize, self->fp))
            break;
    }

    // Tokenize, reallocating if necessary
    int n = strcount(self->buff, sep);
    if(n+1 > self->maxtokens) {
        Tokenizer_free(self->tkz);
        self->maxtokens = n+1;
        self->tkz = Tokenizer_new(self->maxtokens);
        if(self->tkz == NULL)
            return ENOMEM;
    }
    Tokenizer_split(self->tkz, self->buff, sep);
    Tokenizer_strip(self->tkz, extraneous);

    return 0;
}

int GetLineTok_ntokens(GetLineTok *self) {
    return Tokenizer_ntokens(self->tkz);
}

char *GetLineTok_token(Tokenizer * t, int index) {
    return Tokenizer_token(self->tkz, index);
}



    

