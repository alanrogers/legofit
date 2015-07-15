/**
 * @file jobqueue.c
 * @author Alan R. Rogers
 * @brief Multithreaded job queue 
 *
 * This file implements a multithreaded job queue. Jobs are pushed
 * onto a queue by the main program. Each thread (or worker) removes a
 * job from the queue, executes it, and then goes back for
 * another. When all jobs are finished, control returns to the main
 * function. 
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "misc.h"
#include "jobqueue.h"

#if 0
#define DEBUG
#else
#undef DEBUG
#endif

#ifdef DEBUG
#define DPRINTF(arg) printf arg
#else
#define DPRINTF(arg)
#endif

typedef struct Job Job;

/** A single job in the queue */
struct Job {
    /** pointer to next job in queue */
    struct Job *next;

    /** pointer to data for current job */
    void       *param;

    /** pointer to function that undertakes current job */
    int         (*jobfun) (void *param);
};

/** All data used by job queue */
struct JobQueue {
    /** list of jobs */
    Job        *todo;

    /** 0 means not accepting jobs */
    int         acceptingJobs;

    /** number of threads stopped */
    int         threadsStopped;

    /** number of threads */
    int         nthreads;

    /** array of threads */
    pthread_t  *thread;


    /** mutex for locking queue */
    pthread_mutex_t lock;

    /** condition variable for waking workers */
    pthread_cond_t wakeWorker;

    /** condition variable for waking main */
    pthread_cond_t wakeMain;
};

#if 0
pthread_mutex_t stdoutLock = PTHREAD_MUTEX_INITIALIZER;
#endif

void       *threadfun(void *varg);
void        Job_free(Job * job);

JobQueue   *JobQueue_new(int nthreads) {
    int         i, j;
    JobQueue   *jq = malloc(sizeof(JobQueue));

    checkmem(jq, __FILE__, __LINE__);

    jq->todo = NULL;
    jq->acceptingJobs = 1;
    jq->threadsStopped = 0;

    if((i = pthread_mutex_init(&jq->lock, NULL)))
        eprintf("%s:%d: pthread_mutex_init returned %d (%s)",
                __FILE__, __LINE__, i, strerror(i));

    if((i = pthread_cond_init(&jq->wakeWorker, NULL)))
        eprintf("%s:%d: pthread_cond_init returned %d (%s)",
                __FILE__, __LINE__, i, strerror(i));

    if((i = pthread_cond_init(&jq->wakeMain, NULL)))
        eprintf("%s:%d: pthread_cond_init returned %d (%s)",
                __FILE__, __LINE__, i, strerror(i));

    jq->nthreads = nthreads;
    jq->thread = malloc(nthreads * sizeof jq->thread[0]);

    /* launch threads */
    for(i = 0; i < nthreads; ++i) {
        j = pthread_create(&jq->thread[i], NULL, threadfun, jq);
        if(j)
            eprintf("%s:%d: pthread_create returned %d (%s)\n",
                    __FILE__, __LINE__, j, strerror(j));
    }

    return jq;
}

void JobQueue_addJob(JobQueue * jq, int (*jobfun) (void *), void *param) {
    assert(jq);
    int         status;

    Job        *job = malloc(sizeof(Job));

    checkmem(job, __FILE__, __LINE__);

    job->jobfun = jobfun;
    job->param = param;

    status = pthread_mutex_lock(&jq->lock);
    if(status)
        ERR(status, "lock");

    job->next = jq->todo;
    jq->todo = job;

    status = pthread_cond_signal(&jq->wakeWorker);
    if(status)
        ERR(status, "signal wakeWorker");

    status = pthread_mutex_unlock(&jq->lock);
    if(status)
        ERR(status, "unlock");
}

/**
 * Waits until there is a job in the queue, pops it off and executes
 * it, then waits for another.  Runs until jobs are completed and
 * main thread sets acceptingJobs=0.
 */
void       *threadfun(void *vjq) {
    JobQueue   *jq = (JobQueue *) vjq;
    Job        *job;
    int         status;

#ifdef DEBUG
    long unsigned tid = (long unsigned) pthread_self();
#endif

    DPRINTF(("%s %lu entry\n", __func__, tid));

    for(;;) {
        status = pthread_mutex_lock(&jq->lock); /* LOCK */
        if(status)
            ERR(status, "lock");

        while(NULL == jq->todo && 1 == jq->acceptingJobs) {
            DPRINTF(("%s %lu awaiting work. todo=%p acceptingJobs=%d\n",
                     __func__, tid, jq->todo, jq->acceptingJobs));

            status = pthread_cond_wait(&jq->wakeWorker, &jq->lock);
            if(status)
                ERR(status, "wait wakeWorker");
        }

        if(NULL == jq->todo) {
            assert(0 == jq->acceptingJobs);
            DPRINTF(("%s %lu no more jobs\n", __func__, tid));
            break;
        }
        DPRINTF(("%s %lu got work\n", __func__, tid));
        assert(NULL != jq->todo);

        /* remove job from queue */
        job = jq->todo;
        jq->todo = jq->todo->next;

        status = pthread_mutex_unlock(&jq->lock);   /* UNLOCK */
        if(status)
            ERR(status, "unlock");

        DPRINTF(("%s %lu calling jobfun\n", __func__, tid));
        job->jobfun(job->param);
        DPRINTF(("%s %lu back fr jobfun\n", __func__, tid));
        free(job);
    }
    ++jq->threadsStopped;

    if(jq->threadsStopped == jq->nthreads) {
        status = pthread_cond_signal(&jq->wakeMain);
        if(status)
            ERR(status, "signal wakeMain");
    }

    status = pthread_mutex_unlock(&jq->lock);   /* UNLOCK */
    if(status)
        ERR(status, "unlock");

    DPRINTF(("%s %lu exit\n", __func__, tid));
    return NULL;
}

/** Wait until jobs are completed */
void JobQueue_waitOnJobs(JobQueue * jq) {
    int         status;

    DPRINTF(("%s:%d: entry\n", __func__, __LINE__));

    status = pthread_mutex_lock(&jq->lock);
    if(status)
        ERR(status, "lock");

    jq->acceptingJobs = 0;
    while(jq->threadsStopped < jq->nthreads) {
        DPRINTF(("%s:%d: waiting; threadsStopped=%d\n",
                 __func__, __LINE__, jq->threadsStopped));

        status = pthread_cond_broadcast(&jq->wakeWorker);
        if(status)
            ERR(status, "wait wakeWorker");

        status = pthread_cond_wait(&jq->wakeMain, &jq->lock);
        if(status)
            ERR(status, "wait wakeMain");
    }
    DPRINTF(("%s:%d: jobs finished; todo=%p\n",
             __func__, __LINE__, jq->todo));

    status = pthread_mutex_unlock(&jq->lock);
    if(status)
        ERR(status, "unlock");

    DPRINTF(("%s:%d: exit\n", __func__, __LINE__));
}

void Job_free(Job * job) {
    if(NULL == job)
        return;
    Job_free(job->next);
    free(job);
}

void JobQueue_free(JobQueue * jq) {
    assert(jq);
    int         status;

    status = pthread_mutex_destroy(&jq->lock);
    if(status)
        ERR(status, "destroy lock");

    status = pthread_cond_destroy(&jq->wakeWorker);
    if(status)
        ERR(status, "destroy wakeWorker");

    status = pthread_cond_destroy(&jq->wakeMain);
    if(status)
        ERR(status, "destroy wakeMain");

    Job_free(jq->todo);
    free(jq->thread);
    free(jq);
}
