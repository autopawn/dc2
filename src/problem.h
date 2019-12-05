#ifndef DC_PROBLEM_H
#define DC_PROBLEM_H

#include "utils.h"
#include "redstrategy.h"

#define MAX_FILTER 4

#define DEFAULT_TARGET_SOLS 1
#define DEFAULT_THREADS 4
#define DEFAULT_LOCAL_SEARCH 1

typedef enum {
    NO_FILTER = 0,
    BETTER_THAN_EMPTY = 1,
    BETTER_THAN_ONE_PARENT = 2,
    BETTER_THAN_ALL_PARENTS = 3,
    BETTER_THAN_SUBSETS = 4,
} filter;

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

    /* PRECOMPUTATIONS */
    // | Precomputed value of empty solution
    double precomp_empty_value;
    // | Precomputed optimal gain from all clients
    double precomp_client_optimal_gain;
    // | Precomputed distance matrices between facilities (for each mode)
    double **facs_distance[N_FACDIS_MODES];

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

} problem;

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

// Free a problem memory
void problem_free(problem *prob);

// Perform useful precomputations. Must be called before performing operations!
void problem_precompute(problem *prob, redstrategy *rstrats, int n_rstrats);

// Prints a briefing of the problem
void problem_print(const problem *prob, FILE *fp);

#endif