#ifndef DC_LOCALSEARCH_H
#define DC_LOCALSEARCH_H

#include "utils.h"
#include "solution.h"
#include "shuffle.h"

// Performs hill climbing via facility swapings using Witaker's fast exchange heuristic
int solution_whitaker_hill_climbing(const rundata *run, solution *sol, shuffler *shuff);

// Delete repeated solutions on the given array
void solutions_delete_repeated(solution **sols, int *n_sols);

// Perform local searches (in parallel).
void solutions_hill_climbing(rundata *run, solution **sols, int n_sols);

// Updates phi1 and phi2 (arrays with 1st and 2nd nearest solution facility to each client)
// after f_ins is inserted and f_rem is deleted.
void update_phi1_and_phi2(const problem *prob, const solution *sol, int f_ins, int f_rem,
        int *phi1, int *phi2);

#endif