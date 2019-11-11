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

typedef unsigned int uint;

// Auxiliar functions:
void *safe_malloc(size_t size);
uint hash_int(uint x);
void add_to_sorted(int *array, int *len, int val);
void rem_of_sorted(int *array, int *len, int val);

#endif