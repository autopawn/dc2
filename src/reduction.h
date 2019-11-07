#ifndef DC_REDUCTION_H
#define DC_REDUCTION_H

#include "utils.h"
#include "problem.h"
#include "solution.h"

/*
NOTE: All reduction methods expect the solutions to be sorted on decreasing value
*/

typedef enum {
    // | Pick the best solutions
    REDUCTION_BESTS,
    // | Random uniform
    REDUCTION_RANDOM_UNIFORM,
    // | Random with ponderation by rank
    REDUCTION_RANDOM_RANK,
    // | The VR-Heuristic
    REDUCTION_VRHEURISTIC,
    // | Enhanced Glover simple diversity-based starting method
    REDUCTION_GLOVER_SDCE,
    // | EGLOVER_SDCE selecting the best solution of each cluster 
    REDUCTION_GLOVER_SDCE_BESTS,
} reduction_method;

typedef struct {
    const char *nomenclature;
    reduction_method method;
    int n_target;
    int arg;
} redstrategy;

// Allocates an array of redstrategies from nomenclatures 
redstrategy *redstrategy_init_from_nomenclatures(const char **noms, int *n_noms);

// Parses a nomenclature to generate a redstrategy
redstrategy redstrategy_from_nomenclature(const char *nomenclature);

// Apply a redstrategy to reduce a set of solutions
void redstrategy_reduce(const problem *prob, const redstrategy rstrat, 
        solution **sols, int *n_sols);

// Reduce the solutions, just picking the bests
void reduction_bests(const problem *prob, solution **sols, int *n_sols, int n_target);

// Reduce the solutions, picking representatives at random, uniformly.
void reduction_random_uniform(const problem *prob, solution **sols, int *n_sols, int n_target);

// Reduce the solutions, with probabilities according to their rank.
void reduction_random_rank(const problem *prob, solution **sols, int *n_sols, int n_target);


#endif