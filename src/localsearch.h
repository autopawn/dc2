#ifndef DC_LOCALSEARCH_H
#define DC_LOCALSEARCH_H

#include "utils.h"
#include "solution.h"
#include "shuffle.h"


// Fastmat is a big matrix that can be reused over differnt calls of solution_resendewerneck_hill_climbing
typedef struct fastmat fastmat;
fastmat *fastmat_init(int size_y, int size_x);
void fastmat_free(fastmat *mat);

// Performs hill climbing via facility swappings using Resende & Werneck local search
// requires a fastmat of size (prob->n_facs,prob->n_facs) initialized in zeros (retrieves it the same way so it can be reused).
int solution_resendewerneck_hill_climbing(const rundata *run, solution *sol, fastmat *zeroini_mat);

// Performs hill climbing via facility swapings using Witaker's fast exchange heuristic
int solution_whitaker_hill_climbing(const rundata *run, solution *sol, shuffler *shuff);

// Delete repeated solutions on the given array
void solutions_delete_repeated(solution **sols, int *n_sols);

// Perform local searches (in parallel).
void solutions_hill_climbing(rundata *run, solution **sols, int n_sols);

// Updates phi1 and phi2 (arrays with 1st and 2nd nearest solution facility to each client)
// after f_ins is inserted and f_rem is deleted.
void update_phi1_and_phi2(const problem *prob, const solution *sol, int f_ins, int f_rem,
        int *phi1, int *phi2, int *affected);

#endif
