#include "problem.h"

problem *problem_init(int n_facs, int n_clis){
    problem *prob = safe_malloc(sizeof(problem));
    prob->n_facs = n_facs;
    prob->n_clis = n_clis;
    //
    prob->facility_cost = safe_malloc(sizeof(double)*prob->n_facs);
    memset(prob->facility_cost,0,     sizeof(double)*prob->n_facs);
    // Initialize distance cost matrix with rows starting from -1
    prob->distance_cost = safe_malloc(sizeof(double*)*(prob->n_facs+1));
    prob->distance_cost += 1;
    // Initialize row -1
    prob->distance_cost[-1] = safe_malloc(sizeof(double)*prob->n_clis);
    for(int j=0;j<prob->n_clis;j++) prob->distance_cost[-1][j] = INFINITY;
    // Initialize other rrows
    for(int i=0;i<prob->n_facs;i++){
        prob->distance_cost[i] = safe_malloc(sizeof(double)*prob->n_clis);
        memset(prob->distance_cost[i],0,     sizeof(double)*prob->n_clis);
    }

    prob->size_restriction_minimum = -1;
    prob->size_restriction_maximum = -1;

    return prob;
}

problem *problem_copy(const problem *other){
    problem *prob = problem_init(other->n_facs,other->n_clis);
    //
    prob->size_restriction_minimum = other->size_restriction_minimum;
    prob->size_restriction_maximum = other->size_restriction_maximum;
    //
    memcpy(prob->facility_cost,other->facility_cost,sizeof(double)*prob->n_facs);
    //
    for(int i=0;i<prob->n_facs;i++){
        for(int j=0;j<prob->n_clis;j++){
            prob->distance_cost[i][j] = other->distance_cost[i][j];
        }
    }
    return prob;
}

void problem_free(problem *prob){
    // Free facility-client distances
    for(int i=-1;i<prob->n_facs;i++){
        free(prob->distance_cost[i]);
    }
    free(prob->distance_cost-1);
    // Free per facility and client arrays
    free(prob->facility_cost);
    // Free problem
    free(prob);
}

