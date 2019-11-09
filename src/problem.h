#ifndef DC_PROBLEM_H
#define DC_PROBLEM_H

#include "utils.h"

#define N_FACDIS_MODES 3

typedef enum {
    FACDIS_MIN_TRIANGLE = 0,  // df(a,b) = min_j d(a,j)+d(b,j)
    FACDIS_SUM_OF_DELTAS = 1, // df(a,b) = sum_j |d(a,j)-d(b,j)|
    FACDIS_NONE = 2,
} facdismode;

typedef enum {
    SOLDIS_MEAN_SQUARE_ERROR = 0, // D(A,B) = sum_j |d(A,j)-d(B,j)|
    SOLDIS_HAUSDORF = 1,          // D(A,B) = max {sup_a inf_b df(a,b), sup_b inf_a df(b,a)}
    SOLDIS_PER_CLIENT_DELTA = 2,  // D(A,B) = sum_a inf_b df(a,b) + sum_b inf_a df(b,a)
} soldismode;

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

// Perform useful precomputations
void problem_precompute(problem *prob);

// Compute facility-facility distances for a distance mode
void problem_compute_facility_distances(problem *prob, facdismode mode);

// Prints a briefing of the problem
void problem_print(const problem *prob, FILE *fp);

#endif