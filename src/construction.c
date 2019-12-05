#include "construction.h"

// final is expected to have a size >= 2*prob->target_sols 
// cands is liberated.
void update_final_solutions(problem *prob, solution **final, int *n_final,
        solution **cands, int n_cands, int csize){
    // The candidates should have the same number of solutins that csize
    assert(n_cands==0 || cands[0]->n_facs==csize);
    // Perform local search on the candidate solutions
    if(prob->local_search && n_cands>0){
        printf("Performing LS on \033[34;1m%d\033[0m solutions.\n",n_cands);
        solutions_hill_climbing(prob,cands,n_cands);
    }
    // Delete repeated solutions after local search
    int n_cands0 = n_cands;
    solutions_delete_repeated(prob,cands,&n_cands);
    if(n_cands0>n_cands){
        printf("Reduced \033[34;1m%d\033[0m solutions to \033[34;1m%d\033[0m local optima.\n",
            n_cands0,n_cands);
    }
    // Save number of local optima found
    if(prob->local_search){
        prob->lastr_per_size_n_local_optima[csize] = n_cands;
    }
    // Pick the best prob->target_sols candidates
    reduction_bests(prob,cands,&n_cands,prob->target_sols);
    // Merge candidates with the final solutions
    for(int i=0;i<n_cands;i++){
        final[*n_final] = cands[i];
        *n_final += 1;
    }
    // Pick the best prob->target_sols for final solutions
    reduction_bests(prob,final,n_final,prob->target_sols);
    // Update the lower bound
    if(final[0]->value > prob->lower_bound) prob->lower_bound = final[0]->value;
    // Solution cands is not longer needed.
    free(cands);
}

solution **new_find_best_solutions(problem *prob, redstrategy *rstrats, int n_rstrats,
        int *out_n_sols){

    // Set random seed
    srand(prob->random_seed);
    
    // The final solutions:
    int final_n_sols = 0;
    solution **final_sols = safe_malloc(sizeof(solution *)*prob->target_sols*2);
    

    for(int r=0;r<prob->n_restarts;r++){

        printf("\n");
        printf("== RESTART %d/%d ==\n",r+1,prob->n_restarts);

        // The previous generation
        int prev_n_sols = 1;
        solution **prev_sols = safe_malloc(sizeof(solution *)*prev_n_sols);
        prev_sols[0] = solution_empty(prob);

        int csize = 0; // Last base computed
        while(prev_n_sols>0){
            
            printf("\n");
            printf("Base has \033[31;1m%d\033[0m solutions of size \033[32;1m%d\033[0m.\n",prev_n_sols,csize);
            // Apply branch and bound
            if(prob->branch_and_bound){
                int n_sols0 = prev_n_sols;
                branch_and_bound(prob,prev_sols,&prev_n_sols);
                if(prev_n_sols < n_sols0){
                    printf("Pruned \033[31;1m%d\033[0m -> \033[31;1m%d\033[0m solutions, by B&B.\n",n_sols0,prev_n_sols);
                }
            }
            
            // Save number of solutions after expansion
            if(prob->local_search){
                prob->lastr_per_size_n_sols[csize] = prev_n_sols;
            }

            // Apply the reduction strategies
            for(int i=0;i<n_rstrats;i++){
                reduce_by_redstrategy(prob,rstrats[i],prev_sols,&prev_n_sols);
            }


            // Save number of solutions after reduction
            if(prob->local_search){
                prob->lastr_per_size_n_sols_after_red[csize] = prev_n_sols;
            }

            int next_n_sols = 0;
            solution **next_sols = NULL;

            // Expand solutions from the previous generation
            if(csize<prob->n_facs && prev_n_sols>0){
                if(prob->size_restriction_maximum==-1 || csize<prob->size_restriction_maximum){
                    printf("Expanding \033[31;1m%d\033[0m solutions.\n",prev_n_sols);
                    next_n_sols = prev_n_sols;
                    next_sols = new_expand_solutions(prob,prev_sols,prev_n_sols,&next_n_sols);
                }
            }


            // Put the prev generation after LS in the final solutions
            if(prob->size_restriction_minimum==-1 || csize>=prob->size_restriction_minimum){
                update_final_solutions(prob,final_sols,&final_n_sols,prev_sols,prev_n_sols,csize);
            }else{
                for(int i=0;i<prev_n_sols;i++) solution_free(prev_sols[i]);
                free(prev_sols);
            }

            // Now the current gen is the previous one
            prev_n_sols = next_n_sols;
            prev_sols = next_sols;

            // Increase csize
            csize += 1;
            // Update number of iterations
            prob->lastr_n_iterations = csize;
        }

        printf("\n");
    }

    // Retrieve the final solutions:
    *out_n_sols = final_n_sols;
    return final_sols;
}