/// This file maintains an external variable, which is used
/// by all threads to initialize local random number generators.
#include <pthread.h>
#include <gsl/gsl_sf_gamma.h>
pthread_mutex_t seedLock = PTHREAD_MUTEX_INITIALIZER;
unsigned long rngseed=0;
