/**
 * @file xjobqueue.c
 * @author Alan R. Rogers
 * @brief Test jobqueue.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "jobqueue.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

typedef struct {
    double      arg, result;
} TstParam;

typedef struct {
	int i;
} ThreadState;

void *ThreadState_new(void *dat);
void ThreadState_free(void *self);
static void unitTstResult(const char *facility, const char *result);

static void unitTstResult(const char *facility, const char *result) {
    printf("%-26s %s\n", facility, result);
}

void *ThreadState_new(void *dat) {
	ThreadState *ts = malloc(sizeof *ts);
    if(ts==NULL) {
        fprintf(stderr,"%s:%d: bad malloc\n", __FILE__,__LINE__);
        exit(1);
    }

	int *ip = (int *) dat;
	ts->i = *ip;
	return (void *) ts;
}

void ThreadState_free(void *self) {
	free(self);
}

int jobfunc(void *p, void * tdat);

int jobfunc(void *p, void * tdat) {
    TstParam   *param = (TstParam *) p;
	ThreadState *ts = (ThreadState *) tdat;

    param->result = (param->arg) * ts->i;

    return 0;
}

int main(int argc, char **argv) {

    int         verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xjobqueue [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xjobqueue [-v]\n");
        exit(1);
    }

    int         i, njobs = 6, nthreads = 3;
    TstParam    jobs[njobs];
	int         multiplier = 3;
    JobQueue   *jq = JobQueue_new(nthreads,
								  &multiplier,
								  ThreadState_new,
								  ThreadState_free);

    for(i = 0; i < njobs; ++i) {
        jobs[i].arg = i + 1.0;
        jobs[i].result = -99.0;
        JobQueue_addJob(jq, jobfunc, jobs + i);
    }

    JobQueue_waitOnJobs(jq);

    for(i = 0; i < njobs; ++i) {
        if(verbose) {
            printf("%d: %lg --> %lg\n", i, jobs[i].arg, jobs[i].result);
            fflush(stdout);
        }
        assert(jobs[i].result == (i + 1.0) * multiplier);
    }

    for(i = 0; i < njobs; ++i) {
        jobs[i].arg = i + 11.0;
        jobs[i].result = -99.0;
        JobQueue_addJob(jq, jobfunc, jobs + i);
    }

    JobQueue_waitOnJobs(jq);
    JobQueue_noMoreJobs(jq);

    for(i = 0; i < njobs; ++i) {
        if(verbose) {
            printf("%d: %lg --> %lg\n", i, jobs[i].arg, jobs[i].result);
            fflush(stdout);
        }
        assert(jobs[i].result == (i + 11.0) * multiplier);
    }

    JobQueue_free(jq);
    unitTstResult("JobQueue", "OK");
    return 0;
}
