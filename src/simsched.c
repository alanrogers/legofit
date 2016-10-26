#include "simsched.h"
#include "misc.h"
#include <stdio.h>
#include <pthread.h>
#include <string.h>

typedef struct SimSchedLink SimSchedLink;

static SimSchedLink *SimSchedLink_append(SimSchedLink * self, int nOptItr,
                                         int nSimReps);
static SimSchedLink *SimSchedLink_free(SimSchedLink * self);
static SimSchedLink *SimSchedLink_popHead(SimSchedLink * self);
static inline int SimSchedLink_getSimReps(SimSchedLink * self);
static inline int SimSchedLink_getOptItr(SimSchedLink * self);

struct SimSchedLink {
    SimSchedLink *next;
    int         nOptItr, nSimReps;
};

struct SimSched {
    pthread_mutex_t lock;
    SimSchedLink *list;
};

// Append a new link to the end of the list. Return a pointer to the
// beginning of the list.
static SimSchedLink *SimSchedLink_append(SimSchedLink * self, int nOptItr,
                                         int nSimReps) {
    if(self == NULL) {
        SimSchedLink *new = malloc(sizeof(SimSchedLink));
        CHECKMEM(new);
        new->next = NULL;
        new->nOptItr = nOptItr;
        new->nSimReps = nSimReps;
        return new;
    }
    self->next = SimSchedLink_append(self->next, nOptItr, nSimReps);
    return self;
}

// Free the list and return NULL.
static SimSchedLink *SimSchedLink_free(SimSchedLink * self) {
    if(self == NULL)
        return NULL;
    self->next = SimSchedLink_free(self->next);
    free(self);
    return NULL;
}

// Delete head of the list and return a pointer
// to the new head.
static SimSchedLink *SimSchedLink_popHead(SimSchedLink * self) {
    if(self == NULL)
        return NULL;
    SimSchedLink *head = self->next;
    free(self);
    return head;
}

static inline int SimSchedLink_getSimReps(SimSchedLink * self) {
    if(self == NULL) {
        fprintf(stderr, "%s:%d: NULL SimSchedLink\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    return self->nSimReps;
}

static inline int SimSchedLink_getOptItr(SimSchedLink * self) {
    if(self == NULL) {
        fprintf(stderr, "%s:%d: NULL SimSchedLink\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    return self->nOptItr;
}

// Allocate a new SimSched with one stage.
SimSched   *SimSched_new(int nOptItr, int nSimReps) {
    SimSched   *self = malloc(sizeof(SimSched));
    CHECKMEM(self);

    self->list = SimSchedLink_append(NULL, nOptItr, nSimReps);
    CHECKMEM(self->list);

    int         status;

    if((status = pthread_mutex_init(&self->lock, NULL)))
        ERR(status, "lock");

    return self;
}

// Append a stage to a SimSched.
SimSched   *SimSched_append(SimSched * self, int nOptItr, int nSimReps) {
    int         status;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    self->list = SimSchedLink_append(self->list, nOptItr, nSimReps);
    CHECKMEM(self->list);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");

    return self;
}

// Free a SimSched.
void SimSched_free(SimSched * self) {

    self->list = SimSchedLink_free(self->list);
    free(self);
}

int SimSched_getSimReps(SimSched * self) {
    int         status, nSimReps;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    nSimReps = SimSchedLink_getSimReps(self);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");

    return nSimReps;
}

int SimSched_getOptItr(SimSched * self) {
    int         status, nOptItr;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    nOptItr = SimSchedLink_getOptItr(self);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");

    return nOptItr;
}

// If we're at the last stage already, return 1. Otherwise,
// advance to the next stage and return 0.
int SimSched_next(SimSched * self) {
    int         status;

    if(self->list == NULL)
        return 1;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    self->list = SimSchedLink_popHead(self->list);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");

    if(self->list)
        return 1;

    return 0;
}

