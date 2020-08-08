#ifndef DC_FASTMAT_H
#define DC_FASTMAT_H

/* A fastmatrix is a matrix that allows a fast iteration over non-zero elements.
It is build as a dynamic hash table
*/

#include "lib/cadts_hashtable.h"

// coordinate in the fastmat
typedef struct {
    int x,y;
} coord;

// cell in the fastmat
typedef struct {
    // current cell value
    double value;
    // number of clients adding to this value (when 0, the value should be 0)
    int n_clis;
} cell;

CADTS_HASHTABLE(coordcelltable,coord,cell,A.x==B.x && A.y==B.y,A.x+32749*A.y)

typedef struct {
    // Hash table to get index in non-zero array
    coordcelltable *tab;
} fastmat;

// Initialize a fastmat
fastmat *fastmat_init();

// Adds a value on the given position in the fastmatrix
void fastmat_add(fastmat *mat, int y, int x, double v);

// Intended for reverting an addition on the given position in a fastmatrix
void fastmat_rem(fastmat *mat, int y, int x, double v);

// Free fastmat from memory
void fastmat_free(fastmat *mat);

// Resets a fastmat (so that it only contains zeroes)
void fastmat_clean(fastmat *mat);

// Retrieves a value at the given position
double fastmat_get(fastmat *mat, int x, int y);

#endif
