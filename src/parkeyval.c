#include "parkeyval.h"
#include "misc.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#define MAX_PARAM_NAME 100

struct ParKeyVal {
    char        key[MAX_PARAM_NAME];
    double     *valPtr;
    ParKeyVal  *next;
};

ParKeyVal  *ParKeyVal_new(const char *key, double *vptr, ParKeyVal *next);

// Constructor. Call with next=NULL to terminate linked list.
ParKeyVal  *ParKeyVal_new(const char *key, double *vptr, ParKeyVal * next) {
    if(strlen(key) >= MAX_PARAM_NAME)
        eprintf("%s:%s:%d: Parameter name too long. Max=%d.\n",
                __FILE__, __func__, __LINE__, MAX_PARAM_NAME);
    ParKeyVal  *self = malloc(sizeof(ParKeyVal));
    CHECKMEM(self);
    snprintf(self->key, sizeof(self->key), "%s", key);
    self->valPtr = vptr;
    self->next = next;
    return self;
}

// Destructor
void ParKeyVal_free(ParKeyVal * node) {
    if(node == NULL)
        return;
    ParKeyVal_free(node->next);
    free(node);
}

/// Insert a new key/pointer pair into sorted linked list.
ParKeyVal  *ParKeyVal_add(ParKeyVal * node, const char *key, double *vptr) {
    if(node == NULL)
        return ParKeyVal_new(key, vptr, NULL);

    int         i = strcmp(key, node->key);
    if(i < 0)
        return ParKeyVal_new(key, vptr, node);
    else if(i == 0)
        eprintf("%s:%s:%d: Duplicate key \"%s\".\n",
                __FILE__, __func__, __LINE__, key);

    // else..
    node->next = ParKeyVal_add(node->next, key, vptr);
    return node;
}

/// Find key in linked list. On success, set *valPtr to corresponding
/// value and return 0. On failure, set *valPtr to NULL and return 1.
int ParKeyVal_get(ParKeyVal * node, double **valPtr, const char *key) {

    if(node == NULL) {
        // Failed
        *valPtr = NULL;
        return 1;
    }

    int         i = strcmp(key, node->key);
    if(i < 0) {
        // Failed
        *valPtr = NULL;
        return 1;
    } else if(i == 0) {
        // Success
        *valPtr = node->valPtr;
        return 0;
    }
    // else..
    return ParKeyVal_get(node->next, valPtr, key);
}
