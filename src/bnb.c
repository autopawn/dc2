#include "bnb.h"

void branch_and_bound(rundata *run, solution **sols, int *n_sols){
    // Update lower bound
    for(int i=0;i<*n_sols;i++){
        if(sols[i]->value > run->bnb_lower_bound){
            run->bnb_lower_bound = sols[i]->value;
        }
    }
    // Filters solutions with larger upperbound
    int n_sols2 = 0;
    for(int i=0;i<*n_sols;i++){
        double upbound = solution_upper_bound(run,sols[i]);
        if(upbound < run->bnb_lower_bound){
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
