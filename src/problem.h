#ifndef DC_PROBLEM_H
#define DC_PROBLEM_H

#include "utils.h"
#include "redstrategy.h"

#define MAX_FILTER 4

#define DEFAULT_TARGET_SOLS 1
#define DEFAULT_THREADS 4
#define DEFAULT_LOCAL_SEARCH_ONLY_TERMINAL 1
#define DEFAULT_LOCAL_SEARCH_SIZE_CHANGE_MOVEMENTS_ENABLED 1
#define DEFAULT_BRANCHING_FACTOR (-1)

typedef enum {
    NO_FILTER = 0,
    BETTER_THAN_EMPTY = 1,
    BETTER_THAN_ONE_PARENT = 2,
    BETTER_THAN_ALL_PARENTS = 3,
    BETTER_THAN_SUBSETS = 4,
} filter;

#define FILTER_DEFAULT BETTER_THAN_ALL_PARENTS

typedef enum {
    NO_LOCAL_SEARCH        = 0,
    SWAP_BEST_IMPROVEMENT  = 1, // Whitaker's
    SWAP_FIRST_IMPROVEMENT = 2, // Whitaker's
    SWAP_RESENDE_WERNECK   = 3,
} localsearch;

#define DEFAULT_LOCAL_SEARCH SWAP_BEST_IMPROVEMENT

extern const char *filter_names[];

typedef struct {
    // | Number of facilities and clients.
    int n_facs, n_clis;
    // | Weight of each client.
    double *client_weight;
    // | Cost of each facility.
    double *facility_cost;
    // | Distance matrix between facilities and clients.
    double **distance;
    // | Cost of connecting one unit of weight for each unit of distance.
    double transport_cost;
    // | Amount gained for reaching a unit of weight.
    double client_gain;
    // | Cost of not reaching a unit of weight (generaly 0 or infty), expected to be positive.
    double unassigned_cost;
    // | Unless it is -1, the returned solutions must be of that size or larger.
    int size_restriction_minimum;
    // | Unless it is -1, the returned solutions must be of that size or smaller.
    int size_restriction_maximum;
} problem;

typedef struct {
    // The current problem
    problem *prob;
    // | Filter used after expanding the current solutions.
    filter filter;
    // | If B&B is active
    int branch_and_bound;
    // | Lower bound for the optimal cost (for B&B)
    double lower_bound;
    // | Number of solutions requested
    int target_sols;
    // | Number of threads
    int n_threads;
    // | Perform local search?
    int local_search;
    // | Random seed
    int random_seed;
    // | Restarts
    int n_restarts;
    // | If the local search is just done for terminal nodes
    int local_search_only_terminal;
    // | If the local search has the remove movement
    int local_search_rem_movement;
    // | If the local search has the add movement
    int local_search_add_movement;
    /* Maximum number of solutions that will be generated at random from each solution on the pool.
    if -1, then all solutions will be generated*/
    int branching_factor;
    /* PRECOMPUTATIONS */
    // | Precomputed value of empty solution
    double precomp_empty_value;
    // | Precomputed optimal gain from all clients
    double precomp_client_optimal_gain;
    // | Precomputed distance matrices between facilities (for each mode)
    double **facs_distance[N_FACDIS_MODES];
    // | Precomputed facility indexes by proximity for each client for Resende and Werneck's local search
    int **nearly_indexes;

    // | Verbose mode
    int verbose;

    /* OUTPUT INFO */
    // | Total number of iterations performed. Just the first restart when using restarts.
    int firstr_n_iterations;
    // | Number of solutions after expansion. Just the first restart when using restarts.
    int *firstr_per_size_n_sols;
    // | Number of solutions after reduction. Just the first restart when using restarts.
    int *firstr_per_size_n_sols_after_red;
    // | Number of different local optima found. Just the first restart when using restarts.
    int *firstr_per_size_n_local_optima;
    // | Total number of iterations
    int total_n_iterations;
    // | Total number of local searches
    long long int n_local_searches;
    // | Total number of local search movements
    long long int n_local_search_movements;
    // | CPU time performing local search:
    double local_search_seconds;
    // | Time taken on each restart
    double *restart_times;
    // | Values on each restart
    double *restart_values;
} rundata;

// | Retrieves the value of assigning the client c to the facility f
static inline double problem_assig_value(const problem *prob, int f, int c){
    double cweight = prob->client_weight[c];
    if(f==-1){
        return cweight*(-prob->unassigned_cost);
    }else{
        return cweight*(prob->client_gain-prob->transport_cost*prob->distance[f][c]);
    }
}

// Initializes a problem along with all the needed arrays.
problem *problem_init(int n_facs, int n_clis);

// Initializes a problem copying data from another one.
problem *problem_copy(const problem *other);

// Free a problem memory
void problem_free(problem *prob);

// Creates a rundata for the given problem and performs precomputations
rundata *rundata_init(problem *prob, redstrategy *rstrats, int n_rstrats, int n_restarts, int precomp_nearly_indexes);

// Free a rundata
void rundata_free(rundata *run);

// Prints a briefing of the rundata parameters
void rundata_print(const rundata *run, FILE *fp);

#endif
