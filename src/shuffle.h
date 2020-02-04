#ifndef DC_SHUFFLE_H
#define DC_SHUFFLE_H

#include "utils.h"

// Include ranxonshi256 for random number generation
#define RANXOSHI256_EXTERN
#include "lib/ranxoshi256.h"

typedef struct ranxoshi256 ranxoshi256;

typedef struct {
    int len;
    unsigned int *nums;
    unsigned int n_retrieved;
    ranxoshi256 rngen;
} shuffler;

// Initializes a random shuffler generator that retrieves numbers from 0 len-1
shuffler *shuffler_init(int len);

// Gets the next shuffled number
unsigned int shuffler_next(shuffler *shu);

// Shuffles numbers again, already retrieved numbers are restored
void shuffler_reshuffle(shuffler *shu);

// Free shuffler's memory
void shuffler_free(shuffler *shu);


#endif