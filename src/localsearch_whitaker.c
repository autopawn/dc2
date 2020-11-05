#include "localsearch.h"

#define NO_MOVEMENT (-2)

int solution_whitaker_hill_climbing(const rundata *run, solution **solp, const solution *target, shuffler *shuff){
    solution *sol = *solp;
    const problem *prob = run->prob;

    // Is this first improvement?
    int first_improvement = shuff!=NULL;

    // First and Second nearest facility to each client
    int *phi1 = safe_malloc(sizeof(int)*prob->n_clis);
    int *phi2 = safe_malloc(sizeof(int)*prob->n_clis);
    // Auxiliar array for solution_findout
    double *v = safe_malloc(sizeof(double)*prob->n_facs);
    // Initialize arrays
    for(int i=0;i<prob->n_facs;i++){
        v[i] = -INFINITY;
    }
    // For each client, find the second nearest facility in the solution
    for(int i=0;i<prob->n_clis;i++){
        phi1[i] = sol->assigns[i];
        phi2[i] = solution_client_2nd_nearest(prob,sol,i);
    }

    // Available moves
    availmoves *avail = availmoves_init(prob,sol,target);

    // Best solution found so far (if doing path relinking):
    solution *best_sol = NULL;
    if(avail->path_relinking){
        best_sol = solution_copy(prob,sol);
    }

    // Each movement:
    int n_moves = 0;
    while(avail->n_insertions>0 || avail->n_removals>0){
        // Insertion candidate:
        double best_delta = (avail->path_relinking)? -INFINITY : 0;
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
        int k_end = allow_size_decrease? prob->n_facs+1 : prob->n_facs; // also consider not adding a facility if allow_size_decrease
        for(int k=0;k<k_end;k++){
            int f_ins;
            if(k==prob->n_facs){
                f_ins = -1;
            }else{
                f_ins = first_improvement? shuffler_next(shuff) : k;
            }

            // Ignore insertion if the facility is already present in the solution
            if(f_ins>=0 && !avail->avail_inss[f_ins]) continue;
            // Find best facility to remove after inserting f_ins, and profits
            int f_rem;
            double delta_profit, delta_profit_worem;
            solution_findout(prob,sol,f_ins,v,phi2,avail->avail_rems,
                    &f_rem,&delta_profit,&delta_profit_worem);
            // Update best removal and insertion
            int improvement = 0;
            if(delta_profit>best_delta){
                best_delta = delta_profit;
                best_rem = f_rem;
                best_ins = f_ins;
            }
            if(delta_profit_worem>best_delta && allow_size_increase && (f_ins!=-1 || target==NULL)){
                best_delta = delta_profit_worem;
                best_rem = -1;
                best_ins = f_ins;
            }
            if(best_delta>0) improvement = 1;
            // break the loop on first_improvement
            if(first_improvement && improvement) break;
        }

        // No inserting and no removing is equivalent to not performing a movement
        if(best_ins==-1 && best_rem==-1){
            best_ins = NO_MOVEMENT;
            best_rem = NO_MOVEMENT;
            #ifdef DEBUG
                // Assert that the target was reached
                if(target!=NULL) assert(solutionp_facs_cmp(&sol,&target)==0);
            #endif
        }

        // Stop when no movement results in a better solution:
        if(best_ins==NO_MOVEMENT){
            assert(best_rem==NO_MOVEMENT);
            break;
        }

        // Perform swap:
        assert(best_rem!=NO_MOVEMENT);
        double old_value = sol->value;
        if(best_rem!=-1){
            solution_remove(prob,sol,best_rem,phi2,NULL);
        }
        if(best_ins!=-1){
            solution_add(prob,sol,best_ins,NULL);
        }
        availmoves_register_move(avail,best_ins,best_rem);

        // Update best solution so far (if doing path relinking)
        if(avail->path_relinking){
            assert(best_sol);
            if(best_sol->value < sol->value){
                solution_free(best_sol);
                best_sol = solution_copy(prob,sol);
            }
        }

        assert(avail->path_relinking || sol->value>old_value);
        #ifdef DEBUG
            if(isfinite(best_delta)){
                double error_pred = (sol->value-old_value)-best_delta;
                if(error_pred<0) error_pred *= -1;
                assert(error_pred<1e-5 || isnan(error_pred));
            }
        #endif

        // Update phi1 and phi2
        update_phi1_and_phi2(prob,sol,best_ins,best_rem,phi1,phi2,NULL);

        // Count one move:
        n_moves += 1;
    }

    // Free memory
    availmoves_free(avail);
    free(v);
    free(phi2);
    free(phi1);

    // Set sol to best_sol in case we are doing path relinking
    if(best_sol!=NULL){
        *solp = best_sol;
        #ifdef DEBUG
            // Assert that the target was reached
            if(target!=NULL) assert(solutionp_facs_cmp(&sol,&target)==0);
        #endif
        solution_free(sol);
    }
    return n_moves;
}
