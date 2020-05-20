#include "localsearch.h"

#define NO_MOVEMENT (-2)

int solution_whitaker_hill_climbing(const rundata *run, solution *sol, shuffler *shuff){
    const problem *prob = run->prob;

    // Is this first improvement?
    int first_improvement = run->local_search==SWAP_FIRST_IMPROVEMENT;
    assert(shuff!=NULL || !first_improvement);

    // First and Second nearest facility to each client
    int *phi1 = safe_malloc(sizeof(int)*prob->n_clis);
    int *phi2 = safe_malloc(sizeof(int)*prob->n_clis);
    // Array if a facility appears in the solution
    int *used = safe_malloc(sizeof(int)*prob->n_facs);
    // Auxiliar array for solution_findout
    double *v = safe_malloc(sizeof(double)*prob->n_facs);
    // Initialize arrays
    for(int i=0;i<prob->n_facs;i++){
        used[i] = 0;
        v[i] = -INFINITY;
    }
    for(int i=0;i<sol->n_facs;i++){
        used[sol->facs[i]] = 1;
    }
    // For each client, find the second nearest facility in the solution
    for(int i=0;i<prob->n_clis;i++){
        phi1[i] = sol->assigns[i];
        phi2[i] = solution_client_2nd_nearest(prob,sol,i);
    }
    // Each movement:
    int n_moves = 0;
    while(1){
        // Insertion candidate:
        double best_delta = 0;
        int best_rem = NO_MOVEMENT;
        int best_ins = NO_MOVEMENT;

        // Reshuffle shuffler
        if(first_improvement) shuffler_reshuffle(shuff);

        // Check movements that are allowed and not allowed
        int allow_size_decrease = run->local_search_rem_movement && sol->n_facs>1 &&
            (prob->size_restriction_minimum==-1 || sol->n_facs>prob->size_restriction_minimum);
        int allow_size_increase = run->local_search_add_movement &&
            (prob->size_restriction_maximum==-1 || sol->n_facs<prob->size_restriction_maximum);

        // Consider adding a facility while removing the worst facility
        int k_ini = allow_size_decrease? -1 : 0; // also consider not adding a facility if allow_size_decrease
        for(int k=k_ini;k<prob->n_facs;k++){
            int f_ins;
            if(first_improvement && k!=-1){
                f_ins = shuffler_next(shuff);
            }else{
                f_ins = k;
            }

            // Ignore insertion if the facility is already present in the solution
            if(f_ins>=0 && used[f_ins]) continue;
            // Find best facility to remove after inserting f_ins, and profits
            int f_rem;
            double delta_profit, delta_profit_worem;
            solution_findout(prob,sol,f_ins,v,phi2,&f_rem,&delta_profit,&delta_profit_worem);
            // Update best removal and insertion
            int improvement = 0;
            if(delta_profit>best_delta){
                best_delta = delta_profit;
                best_rem = f_rem;
                best_ins = f_ins;
                improvement = 1;
            }
            if(delta_profit_worem>best_delta && allow_size_increase){
                best_delta = delta_profit_worem;
                best_rem = -1;
                best_ins = f_ins;
                improvement = 1;
            }
            // break the loop on first_improvement
            if(first_improvement && improvement) break;
        }
        // No inserting and no removing is equivalent to not performing a movement
        if(best_ins==-1 && best_rem==-1){
            best_ins = NO_MOVEMENT;
            best_rem = NO_MOVEMENT;
        }

        // Stop when no movement results in a better solution:
        if(best_ins==NO_MOVEMENT) break;

        // Perform swap:
        assert(best_rem!=NO_MOVEMENT);
        double old_value = sol->value;
        if(best_rem!=-1){
            solution_remove(prob,sol,best_rem,phi2,NULL);
            used[best_rem] = 0;
        }
        if(best_ins!=-1){
            solution_add(prob,sol,best_ins,NULL);
            used[best_ins] = 1;
        }
        assert(sol->value>old_value);
        #ifdef DEBUG
            double error_pred = (sol->value-old_value)-best_delta;
            if(error_pred<0) error_pred *= -1;
            assert(error_pred<1e-5 || isnan(error_pred));
        #endif

        // Update phi1 and phi2
        update_phi1_and_phi2(prob,sol,best_ins,best_rem,phi1,phi2,NULL);

        // Count one move:
        n_moves += 1;
    }

    // Free memory
    free(v);
    free(used);
    free(phi2);
    free(phi1);
    return n_moves;
}
