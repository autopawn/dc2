#include "bnb.h"

void branch_and_bound(problem *prob, solution **sols, int *n_sols){
    // Update lower bound
    for(int i=0;i<*n_sols;i++){
        if(sols[i]->value > prob->lower_bound){
            prob->lower_bound = sols[i]->value;
        }
    }
    // Filters solutions with larger upperbound
    int n_sols2 = 0;
    for(int i=0;i<*n_sols;i++){
        double upbound = solution_upper_bound(prob,sols[i]);
        if(upbound < prob->lower_bound){
            // Delete solution
            solution_free(sols[i]);
        }else{
            // Solution survives
            sols[n_sols2] = sols[i];
            n_sols2 += 1;
        }
    }
    *n_sols = n_sols2;
}