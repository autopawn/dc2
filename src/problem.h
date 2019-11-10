#ifndef DC_PROBLEM_H
#define DC_PROBLEM_H

#include "utils.h"
#include "redstrategy.h"

typedef enum {
    NO_FILTER = 0,
    BETTER_THAN_ONE_PARENT = 1,
    BETTER_THAN_ALL_PARENTS = 2,
    BETTER_THAN_SUBSETS = 3,
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
    // | Unless it is -1, the returned solutions must be of that size.
    int size_restriction;
    // | Filter used after expanding the current solutions.
    filter filter;
    // | If B&B is active
    int branch_and_bound;
    // | Lower bound for the optimal cost (for B&B)
    double lower_bound;
    // | Amount of solution requested
    int target_sols;
    // | Precomputed optimal gain from all clients
    double precomp_client_optimal_gain;
    // | Distance matrices between facilities (for each mode)
    double **facs_distance[N_FACDIS_MODES];
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