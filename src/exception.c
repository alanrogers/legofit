#include "exception.h"
#include <stdio.h>
#include <stdlib.h>

Exception *Exception_new(void) {
    Exception *self = malloc(sizeof(Exception));
    if(self==NULL) {
        fprintf(stderr, "%s:%s:%d: bad malloc\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    self->next = NULL;
    self->buff[0] = '\0';
    return self;
}

Exception *Exception_push(Exception *root, Exception *new) {
    Exception *e;
    if(root == NULL)
        return new;
    for(e=root; e; e=e->next)
        ;
    e->next = new;
    return root;
}

void Exception_abort(Exception *self) {
    Exception *e;
    fputs("Exception:\n", stderr);
    for(e=self; e; e=e->next) {
        fputs("  ", stderr);
        fputs(e->buff, stderr);
    }
    Exception_free(self);
    exit(EXIT_FAILURE);
}

void Exception_free(Exception *self) {
    if(self==NULL)
        return;
    Exception_free(self->next);
    free(self);
}
