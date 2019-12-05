#ifndef DC_LOCALSEARCH_H
#define DC_LOCALSEARCH_H

#include "utils.h"
#include "solution.h"

// Performs hill climbing via facility swapings using Witaker's fast exchange heuristic 
int solution_whitaker_hill_climbing(const problem *prob, solution *sol);

// Delete repeated solutions on the given array
void solutions_delete_repeated(const problem *prob, solution **sols, int *n_sols);

// Perform local searches (in parallel).
void solutions_hill_climbing(problem *prob, solution **sols, int n_sols);

#endif