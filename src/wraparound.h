#ifndef ARR_WRAPAROUNDF_H
#define ARR_WRAPAROUNDF_H

#include <stdio.h>
#include <stdlib.h>

typedef struct Wraparound Wraparound;

struct  Wraparound {
    unsigned totsize;     // size as allocated
    unsigned cursize;     // number of entries in buffer
    unsigned curpos;      // index where the next push will be stored
    unsigned char buf[1]; // using struct hack
};

Wraparound *Wraparound_new(unsigned totsize);
void     Wraparound_free(Wraparound *self);
unsigned Wraparound_size  (const Wraparound *self);
void     Wraparound_push(Wraparound *self, const unsigned char c);
unsigned Wraparound_pop(Wraparound *self);
void     Wraparound_print(Wraparound *self, FILE *fp);

#endif
