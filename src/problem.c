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
    prob->target_sols = 10;

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
    // Facility distances not yet computed
    for(int mode=0;mode<N_FACDIS_MODES;mode++){
        prob->facs_distance[mode] = NULL;
    }
    //
    return prob;
}

void problem_free(problem *prob){
    // Free facility-client distances
    for(int i=0;i<prob->n_facs;i++){
        free(prob->distance[i]);
    }
    free(prob->distance);
    // Free facility-facility distance matrices, if they are in use:
    for(int mode=0;mode<N_FACDIS_MODES;mode++){
        if(prob->facs_distance[mode]){
            for(int i=0;i<prob->n_facs;i++){
                free(prob->facs_distance[mode][i]);
            }
            free(prob->facs_distance[mode]);
        }
    }
    // Free per facility and client arrays
    free(prob->facility_cost);
    free(prob->client_weight);
    // Free problem
    free(prob);
}

void problem_precompute(problem *prob){
    // Precompute precomp_client_optimal_gain
    prob->precomp_client_optimal_gain = 0;
    for(int k=0;k<prob->n_clis;k++){
        double best_val = problem_assig_value(prob,-1,k);
        for(int i=0;i<prob->n_facs;i++){
            double val = problem_assig_value(prob,i,k);
            if(best_val < val) best_val = val;
        }
        prob->precomp_client_optimal_gain += best_val;
    }
}

void problem_compute_facility_distances(problem *prob, facdismode mode){
    // Allocate distance matrix between facilities 
    assert(mode<N_FACDIS_MODES);
    prob->facs_distance[mode] = safe_malloc(sizeof(double*)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++){
        prob->facs_distance[mode][i] = safe_malloc(sizeof(double)*prob->n_facs);
    }
    // Allocate distance matrix between facilities 
    for(int a=0;a<prob->n_facs;a++){
        for(int b=0;b<prob->n_facs;b++){
            if(mode==FACDIS_SUM_OF_DELTAS){
                double dist = 0;
                for(int j=0;j<prob->n_clis;j++){
                    double delta = prob->distance[a][j]-prob->distance[b][j];
                    if(delta<0) delta = -delta;
                    dist += delta;
                }
                prob->facs_distance[mode][a][b] = dist;
                prob->facs_distance[mode][b][a] = dist;
            }else if(mode==FACDIS_MIN_TRIANGLE){
                double dist = INFINITY;
                for(int j=0;j<prob->n_clis;j++){
                    double dist_sum = prob->distance[a][j]+prob->distance[b][j];
                    if(dist_sum<dist) dist = dist_sum;
                }
                prob->facs_distance[mode][a][b] = dist;
                prob->facs_distance[mode][b][a] = dist;
            }
        }
    }
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
    fprintf(fp,"# BRANCH_AND_BOUND: %d\n",prob->branch_and_bound);
    fprintf(fp,"# LOWER_BOUND: %lf\n",prob->lower_bound);
    fprintf(fp,"# TARGET_SOLS: %d\n",prob->target_sols);
    fprintf(fp,"# PRECOMP_CLIENT_OPTIMAL_GAIN: %lf\n",prob->precomp_client_optimal_gain);
}