#ifndef DC_PROBLEM_H
#define DC_PROBLEM_H

#include "utils.h"
#include "redstrategy.h"

typedef struct {
    // | Number of facilities and clients.
    int n_facs, n_clis;
    // | Cost of each facility.
    double *facility_cost;
    // | Cost matrix between facilities and clients.
    double **distance_cost;
    // | Unless it is -1, the solutions retrieved must be of this size or larger.
    int size_restriction_minimum;
    // | Unless it is -1, the solutions retrieved must be of this size or smaller.
    int size_restriction_maximum;
} problem;

// | Retrieves the cost of assigning the client c to the facility f
static inline double problem_assig_cost(const problem *prob, int f, int c){
    return prob->distance_cost[f][c]; // NOTE that f can be -1
}
// | Retrieves the value (cost*-1) of assigning the client c to the facility f
static inline double problem_assig_value(const problem *prob, int f, int c){
    return -prob->distance_cost[f][c]; // NOTE that f can be -1
}

// Initializes a problem along with all the needed arrays.
problem *problem_init(int n_facs, int n_clis);

// Initializes a problem copying data from another one.
problem *problem_copy(const problem *other);

// Free a problem memory
void problem_free(problem *prob);

#endif
