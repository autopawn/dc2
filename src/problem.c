#include "problem.h"

const char *filter_names[] = {
    "NO_FILTER",
    "BETTER_THAN_ONE_PARENT",
    "BETTER_THAN_ALL_PARENTS",
    "BETTER_THAN_SUBSETS",
};

problem *problem_init(int n_facs, int n_clis){
    problem *prob = safe_malloc(sizeof(problem));
    memset(prob,0,              sizeof(problem));
    prob->n_facs = n_facs;
    prob->n_clis = n_clis;

    //
    prob->client_weight = safe_malloc(sizeof(double)*prob->n_clis);
    memset(prob->client_weight,0,     sizeof(double)*prob->n_clis);
    //
    prob->facility_cost = safe_malloc(sizeof(double)*prob->n_facs);
    memset(prob->facility_cost,0,     sizeof(double)*prob->n_facs);
    // 
    prob->distance = safe_malloc(sizeof(double*)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++){
        prob->distance[i] = safe_malloc(sizeof(double)*prob->n_clis);
        memset(prob->distance[i],0,     sizeof(double)*prob->n_clis);
    }
    //
    prob->facs_distance = safe_malloc(sizeof(double*)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++){
        prob->facs_distance[i] = safe_malloc(sizeof(double)*prob->n_facs);
        memset(prob->facs_distance[i],0,     sizeof(double)*prob->n_facs);
    }
    return prob;
}

void problem_free(problem *prob){
    //
    for(int i=0;i<prob->n_facs;i++){
        free(prob->distance[i]);
    }
    free(prob->distance);
    //
    for(int i=0;i<prob->n_facs;i++){
        free(prob->facs_distance[i]);
    }
    free(prob->facs_distance);
    // 
    free(prob->facility_cost);
    free(prob->client_weight);
    //
    free(prob);
}

void problem_precompute(problem *prob){
    
}

void problem_print(const problem *prob, FILE *fp){
    fprintf(fp,"== PROBLEM ==\n");
    fprintf(fp,"# N_FACILITIES: %d\n",prob->n_facs);
    fprintf(fp,"# N_CLIENTS: %d\n",prob->n_clis);
    fprintf(fp,"# TRANSPORT_COST: %lf\n",prob->transport_cost);
    fprintf(fp,"# CLIENT_GAIN: %lf\n",prob->client_gain);
    fprintf(fp,"# UNASSIGNED_COST: %lf\n",prob->unassigned_cost);
    fprintf(fp,"# SIZE_RESTRICTION: %d\n",prob->size_restriction);
    fprintf(fp,"# FILTER: %s\n",filter_names[prob->filter]);
    fprintf(fp,"# LOWER_BOUND: %lf\n",prob->lower_bound);
    fprintf(fp,"# BRANCH_AND_BOUND: %d\n",prob->branch_and_bound);
}