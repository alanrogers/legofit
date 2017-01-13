/**
 * @file rngseed.c
 * @author Alan R. Rogers
 * @brief Maintain seed of random number generator
 *
 * This file maintains an external variable, which is used
 * by all threads to initialize local random number generators.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <pthread.h>
pthread_mutex_t seedLock = PTHREAD_MUTEX_INITIALIZER;
unsigned long rngseed=0;
