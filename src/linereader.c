#include "linereader.h"
#include "misc.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct LineReader {
    size_t buffsize;
    char *buff;
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

void LineReader_free(LineReader *self) {
    free(self->buff);
    Tokenizer_free(self->tkz);
    free(self);
}

LineReader *LineReader_new(size_t buffsize) {
    if(buffsize < 2)
        DIE("buffsize argument must be >1");
    LineReader *self = malloc(sizeof(LineReader));
    if(self==NULL)
        DIE("bad malloc");
    self->buffsize = buffsize;
    self->buff = malloc(buffsize);
    if(self->buff == NULL)
        DIE("bad malloc");
    self->buff[0] = '\0';
    return self;
};

/// Return pointer to buff on success, NULL on EOF.
/// Abort if allocation fails.
char *LineReader_read(LineReader *self, FILE *fp) {
    self->buff[0] = '\0';
    do {
        if(NULL == fgets(self->buff, self->buffsize, fp))
            return NULL;
    } while(strempty(self->buff));

    while(!strchr(self->buff, '\n') && !feof(self->fp)) {
        assert(strlen(self->buff) == self->buffsize-1);
        // buffer overflow: reallocate
        size_t oldsize = strlen(self->buff);
        self->buffsize *= 2;
        char *oldbuff = self->buff;
        self->buff = realloc(self->buff, self->buffsize);
        CHECKMEM(self->buff);
        // Try to read the rest of the current line. If there is
        // nothing more to read, that's OK, because we've already
        // filled the first half of self->buff.
        if(NULL == fgets(self->buff + oldsize,
                         self->buffsize - oldsize,
                         fp))
            break;
    }

    return self->buff;
}


