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

void LineReader_free(LineReader *self) {
    free(self->buff);
    free(self);
}

LineReader *LineReader_new(size_t buffsize) {
    if(buffsize < 2)
        DIE("buffsize argument must be >1");
    LineReader *self = malloc(sizeof(LineReader));
    CHECKMEM(self);
    self->buffsize = buffsize;
    self->buff = malloc(buffsize);
    CHECKMEM(self->buff);
    self->buff[0] = '\0';
    return self;
};

/// Read next line. Return pointer to buff on success, NULL on EOF.
/// Abort if allocation fails.
char *LineReader_next(LineReader *self, FILE *fp) {
    if(NULL == fgets(self->buff, self->buffsize, fp))
        return NULL;

    // check for buffer overflow; reallocate if necessary
    while(!strchr(self->buff, '\n') && !feof(fp)) {
        // buffer overflow
        size_t oldsize = strlen(self->buff);
        assert(oldsize == self->buffsize-1);
        self->buffsize *= 2;
        self->buff = realloc(self->buff, self->buffsize);
        CHECKMEM(self->buff);
        // Try to read the rest of the current line. If there is
        // nothing more to read, that's OK, because we've already
        // filled the first half of self->buff.
        if(NULL == fgets(self->buff + oldsize,
                         self->buffsize - oldsize,
                         fp)) {
            break;
        }
    }

    return self->buff;
}


