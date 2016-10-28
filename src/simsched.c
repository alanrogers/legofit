#include "simsched.h"
#include "misc.h"
#include <stdio.h>
#include <pthread.h>
#include <string.h>

typedef struct SimSchedLink SimSchedLink;

static SimSchedLink *SimSchedLink_append(SimSchedLink * self, long nOptItr,
                                         long nSimReps);
static SimSchedLink *SimSchedLink_dup(const SimSchedLink *self);
static SimSchedLink *SimSchedLink_free(SimSchedLink * self);
static SimSchedLink *SimSchedLink_popHead(SimSchedLink * self);
static inline long SimSchedLink_getSimReps(SimSchedLink * self);
static inline long SimSchedLink_getOptItr(SimSchedLink * self);
void        SimSchedLink_print(const SimSchedLink *self, FILE *fp);

struct SimSchedLink {
    SimSchedLink *next;
    long          nOptItr, nSimReps;
};

struct SimSched {
    pthread_mutex_t lock;
    SimSchedLink *list;
};

// Append a new link to the end of the list. Return a pointer to the
// beginning of the list.
static SimSchedLink *SimSchedLink_append(SimSchedLink * self, long nOptItr,
                                         long nSimReps) {
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

SimSchedLink *SimSchedLink_dup(const SimSchedLink *self) {
    if(self == NULL)
        return NULL;

    SimSchedLink *new = memdup(self, sizeof(SimSchedLink));
    CHECKMEM(new);

    new->next = SimSchedLink_dup(self->next);
    return new;
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

static inline long SimSchedLink_getSimReps(SimSchedLink * self) {
    if(self == NULL) {
        fprintf(stderr, "%s:%d: NULL SimSchedLink\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    return self->nSimReps;
}

static inline long SimSchedLink_getOptItr(SimSchedLink * self) {
    if(self == NULL) {
        fprintf(stderr, "%s:%d: NULL SimSchedLink\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    return self->nOptItr;
}

void        SimSchedLink_print(const SimSchedLink *self, FILE *fp) {
    if(self==NULL)
        return;
    fprintf(fp, " [nOptItr=%ld, nSimReps=%ld]", self->nOptItr, self->nSimReps);
    SimSchedLink_print(self->next, fp);
}

// Allocate a new SimSched with one stage.
SimSched   *SimSched_new(void) {
    SimSched   *self = malloc(sizeof(SimSched));
    CHECKMEM(self);

    self->list = NULL;

    int         status;
    if((status = pthread_mutex_init(&self->lock, NULL)))
        ERR(status, "init");

    return self;
}

int         SimSched_empty(const SimSched *self) {
    return self->list == NULL;
};

SimSched   *SimSched_dup(const SimSched *self) {
    if(self == NULL)
        return NULL;
    SimSched *new = SimSched_new();
    CHECKMEM(new);

    new->list = SimSchedLink_dup(self->list);
    return new;
}

// Append a stage to a SimSched.
void SimSched_append(SimSched * self, long nOptItr, long nSimReps) {
    int         status;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    self->list = SimSchedLink_append(self->list, nOptItr, nSimReps);
    CHECKMEM(self->list);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");
}

// Free a SimSched.
void SimSched_free(SimSched * self) {

    self->list = SimSchedLink_free(self->list);
    free(self);
}

long SimSched_getSimReps(SimSched * self) {
    int         status;
    long        nSimReps;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    nSimReps = SimSchedLink_getSimReps(self->list);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");

    return nSimReps;
}

long SimSched_getOptItr(SimSched * self) {
    int         status;
    long        nOptItr;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    nOptItr = SimSchedLink_getOptItr(self->list);

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

void        SimSched_print(const SimSched *self, FILE *fp) {
    fprintf(fp, "SimSched:");
    SimSchedLink_print(self->list, fp);
    putc('\n', fp);
}
