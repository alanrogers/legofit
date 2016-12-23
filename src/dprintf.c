/**
 * @file dprintf.c
 * @author Alan R. Rogers
 * @brief Define external mutex for locking output.
 */
#include <pthread.h>
pthread_mutex_t outputLock = PTHREAD_MUTEX_INITIALIZER;
