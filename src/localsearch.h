#ifndef DC_LOCALSEARCH_H
#define DC_LOCALSEARCH_H

#include "utils.h"
#include "solution.h"
#include "shuffle.h"

// An availmoves contains 2 array lists with the available insertion and removal movements
typedef struct {
    // Whether the search is moving to another solution
    int path_relinking;
    // Number of available insertions
    int n_insertions;
    // Array with with the available insertions
    int *insertions;
    // Number of available removals
    int n_removals;
    // Array with the available removals
    int *removals;
    // Directly check if a particular insertion is available
    int *avail_inss;
    // Directly check if a particular removal is available
    int *avail_rems;
    // If a facility is present in the solution
    int *used;
} availmoves;

availmoves *availmoves_init(const problem *prob, const solution *sol, const solution *tgt);
void availmoves_register_move(availmoves *av, int f_ins, int f_rem);
void availmoves_free(availmoves *av);

// Fastmat is a big matrix that can be reused over differnt calls of solution_resendewerneck_hill_climbing
typedef struct fastmat fastmat;
fastmat *fastmat_init(int size_y, int size_x);
void fastmat_free(fastmat *mat);

// Performs hill climbing via facility swappings using Resende & Werneck local search
// requires a fastmat of size (prob->n_facs,prob->n_facs) initialized in zeros (retrieves it the same way so it can be reused).
int solution_resendewerneck_hill_climbing(const rundata *run, solution **solp, const solution *target, fastmat *zeroini_mat);

// Performs hill climbing via facility swapings using Whitaker's fast exchange heuristic
// If a shuffler is provided, 1st improvement is assumed.
int solution_whitaker_hill_climbing(const rundata *run, solution **solp, const solution *target, shuffler *shuff);

// Sort the given array and delete repeated solutions
void solutions_sort_and_delete_repeated(solution **sols, int *n_sols);

// Perform local searches (in parallel).
void solutions_hill_climbing(rundata *run, solution **sols, int n_sols);

// Perform path relinking (in parallel).
void solutions_path_relinking(rundata *run, solution ***sols, int *n_sols);


// Updates phi1 and phi2 (arrays with 1st and 2nd nearest solution facility to each client)
// after f_ins is inserted and f_rem is deleted.
void update_phi1_and_phi2(const problem *prob, const solution *sol, int f_ins, int f_rem,
        int *phi1, int *phi2, int *affected);

#endif
