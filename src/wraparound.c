// Sort of a stack. You can push as many items onto it as you want, but it
// saves only the last totsize items. When you pop items off, they appear
// in an arbitrary order.
#include "wraparound.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

Wraparound *Wraparound_new(unsigned totsize) {
    if(totsize==0u || (totsize & (totsize - 1u))) {
        fprintf(stderr,"%s:%s:%d: totsize=%u; must be power of 2\n",
                __FILE__,__func__,__LINE__, totsize);
        exit(0);
    }

    // Using struct hack.
    size_t size = sizeof(Wraparound) + (totsize-1u)*sizeof(unsigned char);
    Wraparound *self = malloc(size);
    if(self==NULL)
        return NULL;

    self->curpos = self->cursize = 0;
    self->totsize = totsize;
    return self;
}

void Wraparound_free(Wraparound *self) {
    free(self);
}

// Return the number of items in the buffer.
unsigned Wraparound_size(const Wraparound *self) {
    return self->cursize;
}

void Wraparound_push(Wraparound *self, const unsigned char c) {
    self->buf[self->curpos] = c;
    self->curpos = (self->curpos + 1u) & (self->totsize - 1u);
    if(self->cursize != self->totsize)
        ++self->cursize;
}

unsigned Wraparound_pop(Wraparound *self) {
    unsigned c;
    if (self->cursize != 0)  {
      c = self->buf[self->cursize-1u];
      --self->cursize;
      self->curpos = self->cursize;
    } else {
        fprintf(stderr,"%s:%d: can't pop an empty Wraparoundf\n",
                __FILE__, __LINE__);
        exit(1);
    }
    return c;
}

void Wraparound_print(Wraparound *self, FILE *fp) {
    fprintf(stderr,"curpos=%u cursize=%u totsize=%u\n",
            self->curpos, self->cursize, self->totsize);
    int i;
    for(i=0; i < self->cursize; ++i)
        fprintf(stderr," %u", self->buf[i]);
    putc('\n', stderr);
}
