
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
#include <stdbool.h>
#define MAX_PARAM_NAME 100

struct ParKeyVal {
    char key[MAX_PARAM_NAME];
    double *valPtr;             // not locally owned
    ParamStatus pstat;
    ParKeyVal *next;
};

ParKeyVal *ParKeyVal_new(const char *key, double *vptr,
                         ParamStatus pstat, ParKeyVal * next);

/// Constructor. Call with next=NULL to terminate linked list.
ParKeyVal *ParKeyVal_new(const char *key, double *vptr,
                         ParamStatus pstat, ParKeyVal * next) {
    if(strlen(key) >= MAX_PARAM_NAME)
        eprintf("%s:%s:%d: Parameter name too long. Max=%d.\n",
                __FILE__, __func__, __LINE__, MAX_PARAM_NAME);
    ParKeyVal *self = malloc(sizeof(ParKeyVal));
    CHECKMEM(self);
    snprintf(self->key, sizeof(self->key), "%s", key);
    self->valPtr = vptr;
    self->pstat = pstat;
    self->next = next;
    ParKeyVal_sanityCheck(self, __FILE__, __LINE__);
    return self;
}

// Destructor
void ParKeyVal_free(ParKeyVal * self) {
    if(self == NULL)
        return;
    ParKeyVal_free(self->next);
    free(self);
}

/// Insert a new key/pointer pair into sorted linked list.
ParKeyVal *ParKeyVal_add(ParKeyVal * self, const char *key,
                         double *vptr, ParamStatus pstat) {
    if(self == NULL)
        return ParKeyVal_new(key, vptr, pstat, NULL);

    int i = strcmp(key, self->key);
    if(i < 0)
        return ParKeyVal_new(key, vptr, pstat, self);
    else if(i == 0)
        eprintf("%s:%s:%d: Duplicate key \"%s\".\n",
                __FILE__, __func__, __LINE__, key);

    // else..
    self->next = ParKeyVal_add(self->next, key, vptr, pstat);
    return self;
}

/// Find key in linked list. On success, return pointer corresponding
/// to key. On failure, return NULL.
double *ParKeyVal_get(ParKeyVal * self, ParamStatus * pstat, const char *key) {
    if(self == NULL)
        return NULL;

    int i = strcmp(key, self->key);
    if(i < 0)                   // Failed
        return NULL;

    if(i == 0) {                // Success
        *pstat = self->pstat;
        return self->valPtr;
    }

    return ParKeyVal_get(self->next, pstat, key);
}

/// Print a ParKeyVal
void ParKeyVal_print(ParKeyVal * self, FILE * fp) {
    if(self == NULL)
        fprintf(fp, "NULL\n");
    else {
        fprintf(fp, "[%p:%s,%p,", self, self->key, self->valPtr);
        switch (self->pstat) {
        case Free:
            fputs("Free", fp);
            break;
        case Fixed:
            fputs("Fixed", fp);
            break;
        case Gaussian:
            fputs("Gaussian", fp);
            break;
        case Constrained:
            fputs("Constrained", fp);
            break;
        default:
            fprintf(stderr, "%s:%d: Unknown ParamStat value: %d.\n",
                    __FILE__, __LINE__, self->pstat);
            exit(EXIT_FAILURE);
        }
        fputs("]->", fp);
        ParKeyVal_print(self->next, fp);
    }
}

/// Abort if ParKeyVal fails tests
void ParKeyVal_sanityCheck(ParKeyVal * self, const char *file, int line) {
#ifndef NDEBUG
    if(self == NULL)
        return;
    REQUIRE(self->valPtr != NULL, file, line);
    ParKeyVal_sanityCheck(self->next, file, line);
#endif
}

/// Return 1 if the two linked lists are equal; 0 if they differ.
int ParKeyVal_equals(ParKeyVal * lhs, ParKeyVal * rhs) {
    if(lhs == NULL && rhs == NULL)
        return 1;
    if(lhs == NULL && rhs != NULL) {
        return 0;
    }
    if(lhs != NULL && rhs == NULL) {
        return 0;
    }
    assert(lhs != NULL && rhs != NULL);
    if(0 != strcmp(lhs->key, rhs->key)) {
        return 0;
    }
    if(lhs->valPtr == NULL && rhs->valPtr != NULL) {
        return 0;
    }
    if(lhs->valPtr != NULL && rhs->valPtr == NULL) {
        return 0;
    }
    if(lhs->valPtr != NULL && rhs->valPtr != NULL) {
        if(!Dbl_equals_allowNonfinite(*lhs->valPtr, *rhs->valPtr)) {
            return 0;
        }
        if(lhs->pstat != rhs->pstat) {
            return 0;
        }
    }
    return ParKeyVal_equals(lhs->next, rhs->next);
}

/// Return nonzero if name is legal; zero otherwise
int legalName(const char *name) {
    const char *legal =
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ" "0123456789" "._:@$";
    return strlen(name) == strspn(name, legal);
}
