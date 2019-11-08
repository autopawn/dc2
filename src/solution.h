#ifndef DC_SOLUTION_H
#define DC_SOLUTION_H

#include "utils.h"
#include "problem.h"

typedef struct {
    int n_facs;
    // ^ Number of facilities (size) of this solution.
    int *facs;
    // ^ Indexes of the facilities. Sorted.
    int *assigns;
    // ^ For each client, which facility it is assigned to. -1 means unnasigned.
    double value;
    // ^ Value of the solution (> is better)
} solution;

// Solution comparison to sort solution pointers on decreasing value
int solutionp_value_cmp_inv(const void *a, const void *b);

// Creates a new, empty solution.
solution *solution_empty(const problem *prob);

// Creates a solution copying another
solution *solution_copy(const problem *prob, const solution *sol);

// Add a facility to an existing solution
void solution_add(const problem *prob, solution *sol, int newf);

// An upper bound for the best value that a children solution could have
double solution_upper_bound(const problem *prob, const solution *sol);

// Delete solution
void solution_free(solution *sol);

// Compute dissimilitude between solutions based 
double solution_dissimilitude(const problem *prob,
        const solution *sol1, const solution *sol2,
        soldismode sdismode, facdismode fdismode);

// Print a solution to the given descriptor
void solution_print(const problem *prob, const solution *sol, FILE *fp);

#endif