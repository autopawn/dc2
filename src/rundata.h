#ifndef DC2_RUNDATA_H
#define DC2_RUNDATA_H

#include "problem.h"
#include "runinfo.h"
#include "runprecomp.h"

#define MAX_FILTER 4

#define DEFAULT_TARGET_SOLS 1
#define DEFAULT_THREADS 4
#define DEFAULT_LOCAL_SEARCH_SIZE_CHANGE_MOVEMENTS_ENABLED 1
#define DEFAULT_BRANCHING_FACTOR (-1)
#define DEFAULT_BRANCHING_CORRECTION 1
#define BRANCH_AND_BOUND_DEFAULT 0
#define DEFAULT_LOCAL_SEARCH_BEFORE_SELECT 1
#define DEFAULT_SELECT_ONLY_TERMINAL 1

// Possible filters after child solutions are created
typedef enum {
    NO_FILTER = 0,
    BETTER_THAN_EMPTY = 1,
    BETTER_THAN_ONE_PARENT = 2,
    BETTER_THAN_ALL_PARENTS = 3,
    BETTER_THAN_SUBSETS = 4,
} filter;

#define FILTER_DEFAULT BETTER_THAN_ALL_PARENTS

// Local search types
typedef enum {
    NO_LOCAL_SEARCH        = 0,
    SWAP_BEST_IMPROVEMENT  = 1, // Whitaker's
    SWAP_FIRST_IMPROVEMENT = 2, // Whitaker's
    SWAP_RESENDE_WERNECK   = 3,
} localsearch;

#define DEFAULT_LOCAL_SEARCH SWAP_BEST_IMPROVEMENT

// Path relinking use after each iteration of the algorithm
typedef enum {
    NO_PATH_RELINKING              = 0,
    PATH_RELINKING_1_ITER          = 1,
    PATH_RELINKING_UNTIL_NO_BETTER = 2,
} pathrelinkingmode;

#define DEFAULT_PATH_RELINKING NO_PATH_RELINKING



typedef struct {
    // The current problem being solved
    problem *prob;
    // | Filter used after expanding the current solutions.
    filter filter;
    // | If B&B is active
    int branch_and_bound;
    // | Lower bound for the optimal cost (for B&B)
    double bnb_lower_bound;
    // | Number of solutions requested
    int target_sols;
    // | Number of threads
    int n_threads;
    // | Which local search to perform, if any.
    localsearch local_search;
    // | Which local search to use in path relinking
    localsearch local_search_pr;
    // | Random seed
    int random_seed;
    // | Restarts
    int n_restarts;
    // | If only terminal solutions are selected
    int select_only_terminal;
    // | If the local search is perform before final selection
    int local_search_before_select;
    // | If the local search has the remove movement
    int local_search_rem_movement;
    // | If the local search has the add movement
    int local_search_add_movement;
    /* Maximum number of solutions that will be generated at random from each solution on the pool.
    if -1, then all solutions will be generated
    if 0, then ceil(log2(m/p)) solutions will be generated */
    int branching_factor;
    // Whether to increase the branching factor for the first generations because they start with a pool of size 0
    int branching_correction;
    // If PR is enabled
    pathrelinkingmode path_relinking;

    // | Verbose mode
    int verbose;

    // | Run info
    runinfo *run_inf;

    // | Precomputations
    runprecomp *precomp;

} rundata;

// Creates a rundata for the given problem and performs precomputations
// if precomp_nearly_indexes: Precompute, for each client facilities indexes by proximity
rundata *rundata_init(problem *prob, redstrategy *rstrats, int n_rstrats, int n_restarts, int precomp_nearly_indexes, int n_threads, int verbose);

// Free a rundata
void rundata_free(rundata *run);

// Prints a briefing of the rundata parameters
void rundata_print(const rundata *run, FILE *fp);

#endif
