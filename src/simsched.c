/**
 * @file simsched.c
 * @author Alan R. Rogers
 * @brief Simulation schedule.
 *
 * SimSched implements a linked list of stages. Each stage establishes
 * a limit on the number of optimizer iterations and on the number of
 * simulation replicates.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "simsched.h"
#include "misc.h"
#include <stdio.h>
#include <pthread.h>
#include <string.h>

typedef struct Stage Stage;

static Stage *Stage_append(Stage * self, long nOptItr,
                                         long nSimReps);
static Stage *Stage_free(Stage * self);
static Stage *Stage_popHead(Stage * self);
static inline long Stage_getSimReps(Stage * self);
static inline long Stage_getOptItr(Stage * self);

/// A single stage in the optimization process.
struct Stage {
    Stage        *next;     ///< pointer to next Stage in list
    long          nOptItr;  ///< number of optimizer iterations
    long          nSimReps; ///< number of simulation replicates
};

/// Holds a linked list of Stages and a mutex lock.
struct SimSched {
    pthread_mutex_t lock;
    Stage *list;
};

/// Append a new link to the end of the list. Return a pointer to the
/// beginning of the list.
static Stage *Stage_append(Stage * self, long nOptItr,
                                         long nSimReps) {
    if(self == NULL) {
        Stage *new = malloc(sizeof(Stage));
        CHECKMEM(new);
        new->next = NULL;
        new->nOptItr = nOptItr;
        new->nSimReps = nSimReps;
        return new;
    }
    self->next = Stage_append(self->next, nOptItr, nSimReps);
    return self;
}

/// Free the list and return NULL.
static Stage *Stage_free(Stage * self) {
    if(self == NULL)
        return NULL;
    self->next = Stage_free(self->next);
    free(self);
    return NULL;
}

/// Delete head of the list and return a pointer
/// to the new head.
static Stage *Stage_popHead(Stage * self) {
    if(self == NULL)
        return NULL;
    Stage *head = self->next;
    free(self);
    return head;
}

/// Return the number of simulation replicates of current Stage
static inline long Stage_getSimReps(Stage * self) {
    if(self == NULL) {
        fprintf(stderr, "%s:%d: NULL Stage\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    return self->nSimReps;
}

/// Return the number of optimizer iterations of current Stage
static inline long Stage_getOptItr(Stage * self) {
    if(self == NULL) {
        fprintf(stderr, "%s:%d: NULL Stage\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    return self->nOptItr;
}

/// Allocate a new SimSched with no stages.
SimSched   *SimSched_new(void) {
    SimSched   *self = malloc(sizeof(SimSched));
    CHECKMEM(self);

    self->list = NULL;

    int         status;
    if((status = pthread_mutex_init(&self->lock, NULL)))
        ERR(status, "init");

    return self;
}

/// Append a stage to a SimSched.
void SimSched_append(SimSched * self, long nOptItr, long nSimReps) {
    int         status;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    self->list = Stage_append(self->list, nOptItr, nSimReps);
    CHECKMEM(self->list);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");
}

/// Free a SimSched.
void SimSched_free(SimSched * self) {

    self->list = Stage_free(self->list);
    free(self);
}

/// Return number of simulation reps in current stage
long SimSched_getSimReps(SimSched * self) {
    int         status;
    long        nSimReps;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    nSimReps = Stage_getSimReps(self->list);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");

    return nSimReps;
}

/// Return number of optimizer iterations in current stage
long SimSched_getOptItr(SimSched * self) {
    int         status;
    long        nOptItr;

    status = pthread_mutex_lock(&self->lock);
    if(status)
        ERR(status, "lock");

    nOptItr = Stage_getOptItr(self->list);

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

    self->list = Stage_popHead(self->list);

    status = pthread_mutex_unlock(&self->lock);
    if(status)
        ERR(status, "unlock");

    if(self->list)
        return 1;

    return 0;
}

/// Print SimSched
void        SimSched_print(const SimSched *self, FILE *fp) {
    fprintf(fp, "# %5s %7s %8s\n", "Stage", "nOptItr", "nSimReps");
    Stage *stage;
    int i=0;
    for(stage=self->list; stage != NULL; stage=stage->next) {
        fprintf(fp, "# %5d %7ld %8ld\n",
                i, stage->nOptItr, stage->nSimReps);
        ++i;
    }
}

/// Return number of stages
int         SimSched_nStages(const SimSched *self) {
    Stage *stage;
    int nstages=0;
    for(stage=self->list; stage != NULL; stage=stage->next)
        ++nstages;
    return nstages;
}
