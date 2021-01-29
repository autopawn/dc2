#ifndef DC_SOLUTION_H
#define DC_SOLUTION_H

#include "utils.h"
#include "rundata.h"

typedef struct {
    int n_facs;
    // ^ Number of facilities (size) of this solution.
    int *facs;
    // ^ Indexes of the facilities. Sorted.
    int *assigns;
    // ^ For each client, which facility it is assigned to. -1 means unnasigned.
    double value;
    // ^ Value of the solution (> is better)
    int terminal;
    // ^ If the solution is a terminal one (didn't generate better childs).
} solution;

// solution* comparison to sort solution pointers on decreasing value
int solutionp_value_cmp_inv(const void *a, const void *b);

// solution* comparison for equality
int solutionp_facs_cmp(const void *a, const void *b);

// Creates a new, empty solution.
solution *solution_empty(const problem *prob);

// Creates a solution copying another
solution *solution_copy(const problem *prob, const solution *sol);

// Add a facility to an existing solution
void solution_add(const problem *prob, solution *sol, int newf, int *affected);

// Remove a facility to an existing solution
// The phi2 array is optional, if not NULL it should contain the second nearest facility
// for each client.
void solution_remove(const problem *prob, solution *sol, int newf, int *phi2, int *affected);

// An upper bound for the best value that a children solution could have
double solution_upper_bound(const rundata *run, const solution *sol);

// Delete solution
void solution_free(solution *sol);

// Compute dissimilitude between solutions with the given dissimilitude mode and facility distance mode
double solution_dissimilitude(const rundata *run,
        const solution *sol1, const solution *sol2,
        soldismode sdismode, facdismode fdismode);

// Find the index of the second nearest facility to the given client, on the solution
int solution_client_2nd_nearest(const problem *prob, const solution *sol, int cli);

// Find the best option for removal if f_ins is inserted to the solution
// NOTE: v must be intialized with -INFINITY and have size equal to prob->n_facs.
// ^ It is always reset to that state before returning.
// frem_allowed is prob->n_facs array that indicates, for each facility in the solution,
//      if it can be removed. NULL allows all facilities.
// *out_frem is the best facility to remove.
// *out_profit is the profit from swapping.
// *out_profit_worem is the profit from adding the solution without removing
void solution_findout(const problem *prob, const solution *sol, int f_ins, double *v,
        const int *phi2, int *frem_allowed,
        int *out_f_rem, double *out_profit, double *out_profit_worem);

// Print a solution to the given descriptor
void solution_print(const problem *prob, const solution *sol, FILE *fp);

// Checks that has solution is properly computed, the value is right and each client is assigned to its nearest facilty
int solution_check_integrity(const problem *prob, const solution *sol);

#endif
