#ifndef DC_UTILS_H
#define DC_UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <limits.h>

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
void *safe_realloc(void *original, size_t size);

uint hash_int(uint x);
void add_to_sorted(int *array, int *len, int val);
void rem_of_sorted(int *array, int *len, int val);
int elem_in_sorted(int *array, int len, int val);

// Retrieves on how many values both sorted arrays differ; in O(len1+len2) time
int diff_sorted(int *arr1, int len1, int *arr2, int len2);

// Semaphore initialization
sem_t *dc_semaphore_init();
// Semaphore destruction
void dc_semaphore_free(sem_t *sem);

// Memory usage
void get_memory_usage(int* currRealMem, int* peakRealMem, int* currVirtMem, int* peakVirtMem);

static inline void detect_errno(){
    if(errno!=0) fprintf(stderr,"ERROR detected: \"%s\"\n",strerror(errno));
    assert(errno==0);
}

#endif
