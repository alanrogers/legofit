/**
 * @file parkeyval.c
 * @author Alan R. Rogers
 * @brief Linked list associating parameter names with pointers to
 * parameter values.  
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
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

/// Find key in linked list. On success, return pointer corresponding
/// to key. On failure, return NULL.
double *ParKeyVal_get(ParKeyVal * node, const char *key) {

    if(node == NULL)
		return NULL;

    int         i = strcmp(key, node->key);
    if(i < 0)        // Failed
		return NULL;

    if(i == 0)  // Success
		return node->valPtr;

    return ParKeyVal_get(node->next, key);
}

void ParKeyVal_print(ParKeyVal *self, FILE *fp) {
	if(self == NULL)
		fprintf(fp,"NULL\n");
	else {
		fprintf(fp,"[%p:%s,%p]->", self, self->key, self->valPtr);
		ParKeyVal_print(self->next, fp);
	}
}
