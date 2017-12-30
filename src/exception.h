#ifndef EXCEPTION_INCLUDED
#define EXCEPTION_INCLUDED
#define EXCEPTION_SIZE 1000
#include "typedefs.h"

struct Exception {
    struct Exception *next;
    char buff[EXCEPTION_SIZE];
};

Exception *Exception_new(void);
Exception *Exception_push(Exception *root, Exception *new);
void Exception_abort(Exception *self);
void Exception_free(Exception *self);

#endif
