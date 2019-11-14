#ifndef DC_REDUCTION_H
#define DC_REDUCTION_H

#include "utils.h"
#include "redstrategy.h"
#include "problem.h"
#include "solution.h"

// NOTE: not all methods are in reduction.c, there are also methods in reduction_*.c files.

// Perform a redunction indicated by the given redstrategy
void reduce_by_redstrategy(problem *prob, const redstrategy rstrat, 
        solution **sols, int *n_sols);

// Reduce the solutions, just picking the bests
void reduction_bests(const problem *prob, solution **sols, int *n_sols, int n_target);

// Reduce the solutions, picking representatives at random, uniformly.
void reduction_random_uniform(const problem *prob, solution **sols, int *n_sols, int n_target);

// Reduce the solutions, with probabilities according to their rank.
void reduction_random_rank(const problem *prob, solution **sols, int *n_sols, int n_target);

// Reduce the solutions, using Glover's simple diversity-based starting
// The best_of_clusters flag is used to indentify if the best of each cluster is selected,
// instead of the centroid.
void reduction_diversity_starting(const problem *prob, solution **sols, int *n_sols,
        int n_target, soldismode soldis, facdismode facdis, int bests_of_clusters);

void reduction_vr_heuristic(const problem *prob, solution **sols, int *n_sols,
        int n_target, soldismode soldis, facdismode facdis, int vision_range);

#endif