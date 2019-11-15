#ifndef DC_UTILS_H
#define DC_UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>

#include <pthread.h>
#include <semaphore.h>

/* As OS X does not support unnamed semaphores,
named semaphores should be used instead. */
#ifdef __APPLE__
    #define NAMED_SEMAPHORES
#endif

#ifdef NAMED_SEMAPHORES
    #include <fcntl.h>
#endif

typedef unsigned int uint;

// Auxiliar functions:
void *safe_malloc(size_t size);
uint hash_int(uint x);
void add_to_sorted(int *array, int *len, int val);
void rem_of_sorted(int *array, int *len, int val);

// Semaphore initialization
sem_t *dc_semaphore_init();
// Semaphore destruction
void dc_semaphore_free(sem_t *sem);

#endif