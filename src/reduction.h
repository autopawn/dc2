#ifndef DC_REDUCTION_H
#define DC_REDUCTION_H

#include "utils.h"
#include "redstrategy.h"
#include "problem.h"
#include "rundata.h"
#include "solution.h"

// NOTE: not all functions are defined in reduction.c, there are also functions in reduction_*.c files.

// Perform a redunction indicated by the given redstrategy
void reduce_by_redstrategy(const rundata *run, const redstrategy rstrat,
        solution **sols, int *n_sols);

// Reduce the solutions, just picking the bests
void reduction_bests(const rundata *run, solution **sols, int *n_sols, int n_target);

// Reduce the solutions, picking representatives at random, uniformly.
void reduction_random_uniform(const rundata *run, solution **sols, int *n_sols, int n_target, int elitist);

// Reduce the solutions, with probabilities according to their rank.
void reduction_random_rank(const rundata *run, solution **sols, int *n_sols, int n_target, int elitist);

// Reduce the solutions, using Glover's simple diversity-based starting
// The best_of_clusters flag is used to indentify if the best of each cluster is selected,
// instead of the centroid.
void reduction_diversity_starting(const rundata *run, solution **sols, int *n_sols,
        int n_target, soldismode soldis, facdismode facdis, int bests_of_clusters);

void reduction_vr_heuristic(const rundata *run, solution **sols, int *n_sols,
        int n_target, soldismode soldis, facdismode facdis, int vision_range);




// Remove the worst fraction of the given solutions
void reduction_remove_worst(solution **sols, int *n_sols, double factor);



#endif
